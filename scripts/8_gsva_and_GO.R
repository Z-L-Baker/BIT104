#### Gene Set Analysis R SCRIPT ####

#### Set-up ####
rm(list=ls())
#Set working directory
workDir = "/mnt/clusters/admiral/data/c21082179/BIT104/New_Analysis"  # working directory for the R session
setwd(workDir)

## Read in the .RData file from the omnibus_method_results.R script :>
load("Filtered_DESeq2_Results.RData")

# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# BiocManager::install("rrvgo")

library(GSVA) #needed to run GSVA
library("DESeq2") #needed for VST and other DESeq2 objects even though the modelling itself is not done here
# this and tidyr needed for pivoting and such
library(tidyr) 
#library(biomaRt) # needed to map uniprot to entrez
#library(clusterProfiler) # needed to download KEGG pathways cleanly
library(tibble) # makes a few transformations easier akin to dplyr and tidyr
library(mgcv) # needed for fitting GAMs
library(broom) # helps with tidying output of statistical models for further processing
library(ggplot2) # for plots
library(grid) # for copying text from console
library(pheatmap)
library("GO.db")
library(rrvgo)
library(org.Hs.eg.db)
library("dplyr")


gene2go_df = read.table("metadata/merged_output_long.tsv", header = F, sep = "\t")
colnames(gene2go_df) = c("GeneID", "GOID")

# Get all the info required for the GO terms as we have like 20k and need to group them for sets... and also simplify this a lot
go_terms = unique(gene2go_df$GOID)

go_names = AnnotationDbi::select(GO.db,
                                  keys = go_terms,
                                  columns = c("TERM", "ONTOLOGY"),
                                  keytype = "GOID")


go_bp = go_names %>% filter(ONTOLOGY == "BP") # this takes us down to 15K
go_bp = merge(go_bp, gene2go_df, by = "GOID")

# So we will need a background set...

# All gene IDs (background)
gene_universe = unique(gene2go_df$GeneID) # 16k genes, will need to drop the ones from the list which are not expressed AT ALL in the data

# Variance Stabilising Transformation (VST) for GSVA:
vsd = vst(ddsTC, blind = FALSE)
head(assay(vsd), 3)
colData(vsd)

# Drop the genes down to what is expressed in the vst
go2genes = go_bp %>%
  filter(GeneID %in% row.names(vsd)) %>%
  group_by(GOID) %>%
  summarise(genes = list(unique(GeneID)),
            n_genes = n_distinct(GeneID),
            .groups="drop")

go2genes_filtered = go2genes %>% filter(n_genes >= 15 & n_genes <= 500)
length(go2genes_filtered$GOID)  # 6465 GO IDs left
test = unique(unlist(go2genes_filtered$genes)) # 13005 genes left

go_table = go2genes_filtered %>% transmute(go_id = GOID, term = GOID, score = n_genes)

# Calculate similarity of the GO terms
simMatrix = calculateSimMatrix(go_table$go_id, orgdb="org.Hs.eg.db", ont="BP", method="Rel") # "Removed 14 terms that were not found in orgdb for BP"

# need a score by which to pick the best/most interesting GO term to represent the clusters they get reduced to.
# will aim to pik the one with the highest variance across samples.
# extract the list of gene sets
go_sets = setNames(go2genes_filtered$genes, go2genes_filtered$GOID)

# compute mean expression per GO term per sample
term_means = sapply(go_sets, function(gene_vec) {
  genes = intersect(gene_vec, rownames(vsd))   # keep only expressed genes
  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(vsd)))
  } else {
    return(colMeans(assay(vsd)[genes, , drop = FALSE], na.rm = TRUE))
  }
})
# swap the orientaton
term_means = t(term_means)
# calculate the variance across samples
term_var = apply(term_means, 1, var, na.rm = TRUE) 


# down to 293 clustered GO terms :> comparable numbers to KEGG. works for me
reduced = reduceSimMatrix(simMatrix, term_var, threshold = 0.7, orgdb="org.Hs.eg.db")

cluster_map = reduced %>%
  group_by(parent, parentTerm, cluster) %>%
  summarise(
    indiv_ids   = paste(unique(go), collapse = "; "),
    indiv_terms = paste(unique(term), collapse = "; "),
    n_indivs  = n(),
    .groups = "drop"
  )

indiv_to_parent = cluster_map %>%
  tidyr::separate_rows(indiv_ids, sep = "; ") %>%
  dplyr::select(parent, indiv_ids) %>%
  rename(GOID = indiv_ids)

parent_genes = go2genes_filtered %>%
  inner_join(indiv_to_parent, by = "GOID") %>%
  group_by(parent) %>%
  summarise(
    genes = list(unique(unlist(genes))),
    .groups = "drop"
  )

gsets = setNames(parent_genes$genes, parent_genes$parent)

######################################################################################################


# Take the expression set
expr_vst = assay(vsd)
# check all samples are in the metadata
stopifnot(all(colnames(expr_vst) == meta$ID))

genes_in_gsets = unique(unlist(gsets))
length(genes_in_gsets)
# 12989 genes this time, much better than KEGG
head(expr_vst)

expr_mat = as.matrix(expr_vst)

# make a gsva param using gaussian because the VST will have made it suitable
params = gsvaParam(exprData = expr_mat,
                    geneSets = gsets,
                    kcdf = "Gaussian", 
                    maxDiff = TRUE
                    )

# now run GSVA :>
gsva_res = gsva(params)

#check results
head(gsva_res)

# File name is messing up the pivotting of the dataframe of results 
meta_clean = meta %>% dplyr::select(-c(FileName, Azo, Cyp, Imid, Pro))

# make a long format dataframe of results
df = as.data.frame(t(gsva_res)) |>
  rownames_to_column("ID") |>
  left_join(meta_clean, by="ID") |>
  pivot_longer(cols = -c(ID,Treatment,Time,Timetreatment,Rep),
               names_to="GO_Cluster", values_to="Score")

df$Time = as.numeric(as.character(df$Time))
df$Rep = as.factor(df$Rep)
df$Treatment = relevel(df$Treatment, ref = "Control")

fit_gam_one = function(d){
  # model GAMs with treatment-specific smoothing + main treatment effect over time 
  # setting 6 as maximum knots (6 time points - which can be computationally intensive but won't sacrifice accuracy) 
  # and using restricted maximum likelihoods
  gam(Score ~ Treatment + s(Time, by=Treatment, k=6), data=d, method="ML")
}
# fit a model per GOcluster (so each treatment will have a spline per GO_Cluster - 
# later we will both visualise this and compare to a treatment averaged spline per treatment and control per cluster)
# essentially doing a LRT like in the overall deseq2 model but this time on each treatment with the same model 
# but then making each treatment share the spline with the control, and seeing if the spline trajectory is significantly different
paths = unique(df$GO_Cluster)
model_fits  = setNames(vector("list", length(paths)), paths)
for(p in paths){
  model_fits[[p]] = tryCatch(fit_gam_one(filter(df, GO_Cluster==p)))
}

# Need to check the model fits to ensure they are actually good, representative hits and not just cramming data in a model and throwing it at the wall
dir.create("GAM_diagnostics")
dir.create("GAM_diagnostics/Treatments/")
dir.create("GAM_diagnostics/Treatments/GO")

for(pathway in names(model_fits)) {
  fit = model_fits[[pathway]]
  
  if(inherits(fit, "try-error") || is.null(fit)) next  # skip failed models
  
  safe_pathway = gsub("[^A-Za-z0-9]", "_", pathway)
  
  pdf(paste0("GAM_diagnostics/Treatments/GO/", safe_pathway ,".pdf"), width=10, height=7)
  
  # Title page per pathway
  plot.new()
  title(main = paste("Diagnostics for:", pathway), cex.main=1.5)
  
  #Capture the GAM results :>
  txt = capture.output(print(summary(fit)))
  grid.newpage()
  grid.text(paste(txt, collapse = "\n"), x=0.5, y=0.5, just=c("centre"), gp=gpar(fontfamily="mono", fontsize=9))
  
  # capture the console text details of the gam.check output
  txt = capture.output(gam.check(fit, rep=1))
  grid.newpage()
  grid.text(paste(txt, collapse = "\n"), x=0.5, y=0.5, just=c("centre"), gp=gpar(fontfamily="mono"))
  
  # plot of smooth terms with the confidence intervals
  plot(fit, rug=TRUE, main=paste("Smooth terms for:", pathway))
  
  dev.off()
}

# double check that none literally failed or anything
for(i in seq_along(model_fits)) {
  fit = model_fits[[i]]
  if(inherits(fit, "try-error") || is.null(fit)){
    print(paste0(fit, "in", i, "is NULL"))
  }
}

# Check if any of them CLEARLY violate issues with the p values of the smoothed terms or high deviations in the residuals 
flags = data.frame(GO_Cluster=names(model_fits), PoorFit=NA)
for(i in seq_along(model_fits)) {
  fit = model_fits[[i]]
  if(inherits(fit, "try-error") || is.null(fit)) next
  k_table = mgcv:::k.check(fit)
  res = residuals(fit, type="deviance")
  
  flags$PoorFit[i] = any(k_table[,"p-value"] < 0.05, abs(res) > 3)
}
bad_fits = subset(flags, flags$PoorFit == T)

# There are 24 noted bad fits so I will check each of these manually.
# OK manually checked all the plots look fine, maybe SLIGHTLY out occassionally but nothing enough to suggest any real issues.
#write.csv(bad_fits, "GAM_diagnostics/poor_fits.csv", row.names=FALSE)

# compute time-change significance from fitted GAMs (aka if the GSVA score changes signifcantly)
time_sig_res = lapply(names(model_fits), function(pathway) {
  fit = model_fits[[pathway]]
  if(inherits(fit, "try-error") || is.null(fit)) return(NULL)
  
  s_summ = summary(fit)$s.table  # edf, Ref.df, F, p-value
  
  tibble(
    GO_Cluster = pathway,
    Treatment = gsub("s\\(Time\\):Treatment", "", rownames(s_summ)),  # remove the prefix
    p_value_time = s_summ[, "p-value"],
    time_sig_flag = ifelse(s_summ[, "p-value"] < 0.05, "T", "F")
  )
})
time_sig_res = bind_rows(time_sig_res)

#########
# function to compare each treatment against control and the treatment for a given cluster (nested model to perform LRT on)
compare_to_control = function(d, control="Control"){
  treatments = setdiff(unique(d$Treatment), control)
  out = list()
  
  # full model
  full = gam(Score ~ Treatment + s(Time, by=Treatment, k=6), data=d, method="ML")
  
  for(trt in treatments){
    # reduced model: merge treatment with control (i.e. force them to share the spline)
    d$Treatment_reduced = d$Treatment
    d$Treatment_reduced[d$Treatment==trt] = control
    
    # run reduced model
    reduced = gam(Score ~ Treatment_reduced + s(Time, by=Treatment_reduced, k=6),
                   data=d, method="ML")
    
    # these are basically the same models so I'm fairly happy they will work but I will check for high deviations in the residuals from the forced spline sharing
    if(max(abs(residuals(reduced))) > 3) {
      warning(paste("High residuals in", trt))
    }
    
    # run an anova with Chi sqaured as an LRT proxy to see if the model is significantly improved by the treatment being modelled separately from the control
    # aka is it significantly different from the control :>
    pval = anova(reduced, full, test="Chisq")$`Pr(>Chi)`[2]
    
    # output the treatment and p value as a tibble
    out[[trt]] = tibble(Treatment=trt, p_value=pval)
  }
  
  bind_rows(out)
}

# run all f the control comparison tests pairwise across all treatments with the control within each cluster (this takes a good min)
pairwise_res = lapply(paths, function(p){
  d = filter(df, GO_Cluster==p)
  tmp = compare_to_control(d, control="Control") |>
    mutate(GO_Cluster=p)
})

# bind the rows, do FDR mitigation and 
# make a flag for if the results for the treatment that were significantly improved from being modelled separately to the control or not which will be helpful for plotting
pairwise_res = bind_rows(pairwise_res) |>
  mutate(p_adj = p.adjust(p_value, "BH"),
         sig_flag = ifelse(p_adj < 0.05, "T", "F"))

# add the significance flags into the gsva dataframe
df_with_sig = df %>%
  left_join(pairwise_res %>% dplyr::select(GO_Cluster, Treatment, sig_flag),
            by = c("GO_Cluster", "Treatment"))

## aslo add the sig flags for trajectory changing over time to the plot :>
df_with_sig = df_with_sig %>%
  left_join(time_sig_res %>% dplyr::select(GO_Cluster, Treatment, time_sig_flag),
            by = c("GO_Cluster", "Treatment"))

################# PLOTTING ##############################

df_with_sig = df_with_sig %>%
  mutate(Treatment = factor(Treatment, levels = c("Control", setdiff(unique(Treatment), "Control"))))

# Define treatment colors (adjust as needed)
treatment_colors = c(
  "Control" = "black", "Azo-Pro" = "#AE4D66", "Azoxy" = "#FF9933",
  "Cyp" = "#00CC33", "Cyp-Pro" = "#2E6666", "Imd-Pro" = "#4866CC",
  "Imid" = "#33CCFF", "Pro" = "#5C0099"
)


df_with_sig = df_with_sig %>%
  mutate(
    FacetGroup = case_when(
      Treatment %in% c("Control", "Pro", "Azoxy", "Azo-Pro") ~ "Azo-Pro Set",
      Treatment %in% c("Control", "Pro", "Cyp", "Cyp-Pro") ~ "Cyp-Pro Set",
      Treatment %in% c("Control", "Pro", "Imid", "Imd-Pro") ~ "Imd-Pro Set",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(FacetGroup))

dir.create("figures/GSVA_Treatment_Splines")
dir.create("figures/GSVA_Treatment_Splines/GO_Clusters")

# Define treatment sets
treatment_sets = list(
  "Azo-Pro Set" = c("Control", "Pro", "Azoxy", "Azo-Pro"),
  "Cyp-Pro Set" = c("Control", "Pro", "Cyp", "Cyp-Pro"),
  "Imd-Pro Set" = c("Control", "Pro", "Imid", "Imd-Pro")
)

for (pathway in unique(df_with_sig$GO_Cluster)) {
  # Subset data for this pathway
  d = df_with_sig %>% filter(GO_Cluster == pathway)
  fit = model_fits[[pathway]]
  if (inherits(fit, "try-error") || is.null(fit)) next
  
  safe_pathway = gsub("[^A-Za-z0-9]", "_", pathway)
  
  # Extract significance flags directly from df_with_sig
  sig_flags = d %>% dplyr::select(Treatment, sig_flag, time_sig_flag) %>% distinct()
  
  # Ensure Control is always solid + thick
  sig_flags = sig_flags %>% mutate(sig_flag = ifelse(Treatment == "Control", "T", sig_flag))
  
  # Treatments significant for points
  signif_treatments = sig_flags %>% filter(sig_flag == "T") %>% pull(Treatment)
  keep_points = unique(c("Control", signif_treatments))
  
  # Create prediction grid for all treatment sets
  pred_grid = do.call(rbind, lapply(names(treatment_sets), function(set_name) {
    expand.grid(
      Time = seq(min(d$Time), max(d$Time), length.out = 100),
      Treatment = treatment_sets[[set_name]]
    ) %>% mutate(FacetGroup = set_name)
  }))
  pred_grid = pred_grid %>% mutate(Treatment = factor(Treatment, levels = levels(d$Treatment)))
  
  # Predict GAM with 95% CIs
  pred = predict(fit, newdata = pred_grid, se.fit = TRUE)
  pred_df = pred_grid %>%
    mutate(
      fit = pred$fit,
      se = pred$se.fit,
      lower = fit - 1.96 * se,
      upper = fit + 1.96 * se
    ) %>%
    # Join significance flags
    left_join(sig_flags, by = "Treatment")
  
  # Asterisks for time-significant treatments
  reference_x = 24
  asterisk_x = 25
  
  asterisk_df = pred_df %>%
    filter(time_sig_flag == "T")
  
  if (nrow(asterisk_df) > 0) {
    asterisk_df = asterisk_df %>%
      group_by(Treatment, FacetGroup, Time) %>% 
      summarise(fit = mean(fit), .groups = "drop") %>% 
      group_by(Treatment, FacetGroup) %>%
      summarise(
        fit_ref = approx(Time, fit, xout = reference_x, rule = 2)$y,  # rule=2 for extrapolation
        .groups = "drop"
      ) %>%
      arrange(FacetGroup, fit_ref) %>%  # sort for stacking
      group_by(FacetGroup) %>%
      mutate(
        # stagger vertically: small offset per treatment
        fit = fit_ref + 0.02 * row_number(),  
        Time = asterisk_x
      )
  } else {
    asterisk_df = tibble(Treatment = character(), FacetGroup = character(),
                         fit = numeric(), Time = numeric())
  }
  
  # Plot
  p = ggplot() +
    geom_ribbon(
      data = pred_df,
      aes(x = Time, ymin = lower, ymax = upper, fill = Treatment),
      alpha = 0.15
    ) +
    geom_line(
      data = pred_df,
      aes(x = Time, y = fit, color = Treatment, 
          linetype = sig_flag, linewidth = sig_flag)
    ) +
    geom_text(
      data = asterisk_df,
      aes(x = Time, y = fit, label = "*", color = Treatment),
      inherit.aes = FALSE,
      size = 8,          # bigger asterisks
      show.legend = FALSE
    ) +
    # Optional: show raw points for significant treatments
    #    geom_point(
    #      data = d %>% filter(Treatment %in% keep_points),
    #      aes(x = Time, y = Score, color = Treatment),
    #      alpha = 0.6, size = 1
    #    ) +
    scale_color_manual(values = treatment_colors) +
    scale_fill_manual(values = treatment_colors, guide = "none") +
    scale_linetype_manual(values = c("T" = "solid", "F" = "dotted"), name = "Diff from Control?") +
    scale_linewidth_manual(values = c("T" = 1, "F" = 0.6), guide = "none") +
    theme_classic(base_size = 12) +
    labs(title = pathway, x = "Time (hours)", y = "GSVA Score", color = "Treatment") +
    facet_wrap(~FacetGroup, ncol = 2, scales = "fixed")
  
  # Save
  ggsave(paste0("figures/GSVA_Treatment_Splines/GO_Clusters/", safe_pathway,".png"), p,
         width = 12, height = 7)
}



#### Check for synergy with GSVA results ####




# Make the interaction binary presence indicators again
df$Azo = factor(ifelse(grepl("^Azo", df$Treatment), 1, 0))
df$Cyp = factor(ifelse(grepl("^Cyp", df$Treatment), 1, 0))
df$Imid = factor(ifelse(grepl("^Im", df$Treatment), 1, 0))
df$Pro = factor(ifelse(grepl("Pro$", df$Treatment), 1, 0))
df$Azo_Pro = with(df, ifelse(Azo==1 & Pro==1, 1, 0))
df$Cyp_Pro = with(df, ifelse(Cyp==1 & Pro==1, 1, 0))
df$Imid_Pro = with(df, ifelse(Imid==1 & Pro==1, 1, 0))

fit_int_gam = function(d){
  # model GAMs of the interactions effects over time 
  # setting 6 as maximum knots (6 time points - which can be computationally intensive but won't sacrifice accuracy) 
  # and using restricted maximum likelihoods
  gam(Score ~
        s(Time, k=6) + # Basically no drugs present - control spline
        s(Time, by=Azo, k=6) +
        s(Time, by=Cyp, k=6) +
        s(Time, by=Imid, k=6) +
        s(Time, by=Pro, k=6) +
        s(Time, by=Azo_Pro, k=6) +
        s(Time, by=Cyp_Pro, k=6) +
        s(Time, by=Imid_Pro, k=6),
      data=d,
      method="ML")
}

# fit a model per GOcluster so each predictor and interaction term will have a spline per GO_Cluster - 
# later we will visualise this, although mostly just the interaction splines (above zero shows synergy and below zero shows antagonism)
# will also do essentially an LRT with a reduced model with no interaction terms again, to see if it shows significant model improvement (aka presence of non additive interactions)
paths = unique(df$GO_Cluster)
int_model_fits  = setNames(vector("list", length(paths)), paths)
for(p in paths){
  int_model_fits[[p]] = tryCatch(fit_int_gam(filter(df, GO_Cluster==p)))
}

# Need to check the model fits to ensure they are actually good, representative hits and not just cramming data in a model and throwing it at the wall
dir.create("GAM_diagnostics")
dir.create("GAM_diagnostics/Interactions/")
dir.create("GAM_diagnostics/Interactions/GO")

for(pathway in names(int_model_fits)) {
  int_fit = int_model_fits[[pathway]]
  
  if(inherits(int_fit, "try-error") || is.null(int_fit)) next  # skip failed models
  
  safe_pathway = gsub("[^A-Za-z0-9]", "_", pathway)
  
  pdf(paste0("GAM_diagnostics/Interactions/GO/", safe_pathway ,".pdf"), width=10, height=7)
  
  # Title page per pathway
  plot.new()
  title(main = paste("Diagnostics for:", pathway), cex.main=1.5)
  
  #Capture the GAM results :>
  txt = capture.output(print(summary(int_fit)))
  grid.newpage()
  grid.text(paste(txt, collapse = "\n"), x=0.5, y=0.5, just=c("centre"), gp=gpar(fontfamily="mono", fontsize=9))
  
  # capture the console text details of the gam.check output
  txt = capture.output(gam.check(int_fit, rep=1))
  grid.newpage()
  grid.text(paste(txt, collapse = "\n"), x=0.5, y=0.5, just=c("centre"), gp=gpar(fontfamily="mono"))
  
  # plot of smooth terms with the confidence intervals
  plot(int_fit, rug=TRUE, main=paste("Smooth terms for:", pathway))
  
  dev.off()
}

# double check that none literally failed or anything
for(i in seq_along(int_model_fits)) {
  int_fit = int_model_fits[[i]]
  if(inherits(int_fit, "try-error") || is.null(int_fit)){
    print(paste0(int_fit, "in", i, "is NULL"))
  }
}

# Check if any of them CLEARLY violate issues with the p values of the smoothed terms or high deviations in the residuals 
int_flags = data.frame(GO_Cluster=names(int_model_fits), PoorFit=NA)
for(i in seq_along(int_model_fits)) {
  int_fit = int_model_fits[[i]]
  if(inherits(int_fit, "try-error") || is.null(int_fit)) next
  int_k_table = mgcv:::k.check(int_fit)
  int_res = residuals(int_fit, type="deviance")
  
  int_flags$PoorFit[i] = any(int_k_table[,"p-value"] < 0.01, any(abs(res)) > 3)
}
int_bad_fits = subset(int_flags, int_flags$PoorFit == T)

### About half the models flag up on this, although not overly concerned as there's likely more noise 
# as usual will check any and all gam diagnostics pdfs before trusting the GAM model results

reduced_interactions = function(d){
  interactions = c("Azo_Pro", "Cyp_Pro", "Imd_Pro")
  out = list()
  
  # full model
  full = gam(Score ~
                s(Time, k=6) + # Basically no drugs present - control spline
                s(Time, by=Azo, k=6) +
                s(Time, by=Cyp, k=6) +
                s(Time, by=Imid, k=6) +
                s(Time, by=Pro, k=6) +
                s(Time, by=Azo_Pro, k=6) +
                s(Time, by=Cyp_Pro, k=6) +
                s(Time, by=Imid_Pro, k=6),
              data=d,
              method="ML")
  
  for(inter in interactions){
    
    # choose which interaction to drop
    keep_terms = c("s(Time,k=6)",
                    "s(Time,by=Azo,k=6)",
                    "s(Time,by=Cyp,k=6)",
                    "s(Time,by=Imid,k=6)",
                    "s(Time,by=Pro,k=6)")
    
    if(inter != "Azo_Pro") keep_terms = c(keep_terms, "s(Time, by=Azo_Pro, k=6)")
    if(inter != "Cyp_Pro") keep_terms = c(keep_terms, "s(Time, by=Cyp_Pro, k=6)")
    if(inter != "Imd_Pro") keep_terms = c(keep_terms, "s(Time, by=Imid_Pro, k=6)")
    
    # build formula
    reduced_formula = as.formula(
      paste("Score ~", paste(keep_terms, collapse=" + "))
    )
    
    reduced = gam(reduced_formula, data=d, method="ML")
    
    # test reduced vs full
    pval = anova(reduced, full, test="Chisq")$`Pr(>Chi)`[2]
    
    out[[inter]] = tibble(Interaction=inter, p_value=pval)
  }
  
  bind_rows(out)
}

# run all of the full model vs reduced model tests
int_pairwise_res = lapply(paths, function(p){
  int_d = filter(df, GO_Cluster==p)
  int_tmp = reduced_interactions(int_d) |>
    mutate(GO_Cluster=p)
})

# bind the rows, do FDR mitigation and 
# make a flag for if the results for the interaction terms which significantly improved the model or not which will be helpful for plotting
int_pairwise_res = bind_rows(int_pairwise_res) |>
  mutate(p_adj = p.adjust(p_value, "BH"),
         sig_flag = ifelse(p_adj < 0.05, "T", "F"))

# Create an interaction pathway
df = df %>%
  mutate(
    Interaction = case_when(
      Azo_Pro == 1 ~ "Azo_Pro",
      Cyp_Pro == 1 ~ "Cyp_Pro",
      Imid_Pro == 1 ~ "Imd_Pro",
      TRUE ~ NA_character_
    )
  )
# add the significance flags into the gsva dataframe
int_df_with_sig = df %>%
  left_join(int_pairwise_res %>% dplyr::select(GO_Cluster, Interaction, sig_flag),
            by = c("GO_Cluster", "Interaction"))


#### Plotting Interactions :> ####

dir.create("figures/GSVA_Interaction_Splines")
dir.create("figures/GSVA_Interaction_Splines/GO_Clusters")

interaction_colors = c(
  "Azo_Pro" = "#AE4D66", "Cyp_Pro" = "#2E6666", "Imd_Pro" = "#4866CC"
)

int_df_with_sig = subset(int_df_with_sig, is.na(int_df_with_sig$Interaction) == F)
int_df_with_sig$Interaction = as.factor(int_df_with_sig$Interaction)
interactions = unique(int_df_with_sig$Interaction)

interaction_map = list(
  "Azo_Pro" = data.frame(
    Pro = factor("0", levels = c("0","1")),
    Azo = factor("0", levels = c("0","1")),
    Cyp = factor("0", levels = c("0","1")),
    Imid = factor("0", levels = c("0","1")),
    Azo_Pro = 1, Cyp_Pro = 0, Imid_Pro = 0
  ),
  "Cyp_Pro" = data.frame(
    Pro = factor("0", levels = c("0","1")),
    Azo = factor("0", levels = c("0","1")),
    Cyp = factor("0", levels = c("0","1")),
    Imid = factor("0", levels = c("0","1")),
    Azo_Pro = 0, Cyp_Pro = 1, Imid_Pro = 0
  ),
  "Imd_Pro" = data.frame(
    Pro = factor("0", levels = c("0","1")),
    Azo = factor("0", levels = c("0","1")),
    Cyp = factor("0", levels = c("0","1")),
    Imid = factor("0", levels = c("0","1")),
    Azo_Pro = 0, Cyp_Pro = 0, Imid_Pro = 1
  )
)

for (pathway in unique(int_df_with_sig$GO_Cluster)) {
  # Subset data for this pathway
  d = int_df_with_sig %>% filter(GO_Cluster == pathway)
  fit = int_model_fits[[pathway]]
  if (inherits(fit, "try-error") || is.null(fit)) next  # skip failed models
  
  safe_pathway = gsub("[^A-Za-z0-9]", "_", pathway)
  
  # Grab p-values for this pathway
  sig_res = int_pairwise_res %>% bind_rows() %>%
    filter(GO_Cluster == pathway) %>%
    mutate(sig_flag = ifelse(p_value < 0.05, "T", "F"))
  
  # Skip if no results found or all are NA
  if (nrow(sig_res) == 0 || all(is.na(sig_res$sig_flag))) next
  
  pred_grid = do.call(rbind, lapply(interactions, function(int_name) {
    base = expand.grid(Time = seq(min(d$Time), max(d$Time), length.out = 100))
    base$Interaction = int_name
    mapping = interaction_map[[int_name]][rep(1, nrow(base)), ]
    cbind(base, mapping)
  }))
  
  pred_grid = pred_grid %>%
    mutate(Interaction = factor(Interaction, levels = levels(d$Interaction)))
  
  # Predict GAM with 95% CIs
  pred = predict(fit, newdata = pred_grid, se.fit = TRUE)
  pred_df = pred_grid %>%
    mutate(
      fit = pred$fit,
      se = pred$se.fit,
      lower = fit - 1.96 * se,
      upper = fit + 1.96 * se
    )
  
  # Add significance flags (per interaction)
  sig_flags = sig_res %>% dplyr::select(Interaction, sig_flag) %>% distinct()
  pred_df = pred_df %>% left_join(sig_flags, by = "Interaction")
  
  # Create plot
  p = ggplot() +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    geom_ribbon(
      data = pred_df,
      aes(x = Time, ymin = lower, ymax = upper, fill = Interaction),
      alpha = 0.15
    ) +
    geom_line(
      data = pred_df,
      aes(x = Time, y = fit, color = Interaction, linetype = sig_flag, linewidth = sig_flag)
    ) +
    scale_color_manual(values = interaction_colors) +
    scale_fill_manual(values = interaction_colors, guide = "none") +
    scale_linetype_manual(
      values = c("T" = "solid", "F" = "dotted"),
      labels = c("T" = "Yes", "F" = "No"),
      name = "Significant"
    ) +
    scale_linewidth_manual(values = c("T" = 1, "F" = 0.6), guide = "none") +
    theme_classic(base_size = 10) +
    labs(
      title = paste0(pathway, " Interaction Effect Dynamics") ,
      x = "Time (hours)",
      y = "GSVA Score",
      color = "Interaction"
    )
  
  ggsave(paste0("figures/GSVA_Interaction_Splines/GO_Clusters/", safe_pathway,".png"), p, width = 9, height = 6)
}


#### Save gsets for comparison with kegg results ####

dir.create("sets")

go_sets = gsets

saveRDS(go_sets, file = "sets/go_sets.rds")

# I also want to save the cluster map for future reference :>
write.table(cluster_map, "sets/go_cluster_map.tsv", row.names= F, sep ="\t", quote = F)

# Save .RData file
save.image(file=paste0("GO_GSVA_Results", ".RData"))
