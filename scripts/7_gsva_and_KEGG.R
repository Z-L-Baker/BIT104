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
#BiocManager::install("GSVA")

library(GSVA) #needed to run GSVA
library("DESeq2") #needed for VST and other DESeq2 objects even though the modelling itself is not done here
library("dplyr") # this and tidyr needed for pivoting and such
library(tidyr) 
library(biomaRt) # needed to map uniprot to entrez
library(clusterProfiler) # needed to download KEGG pathways cleanly
library(tibble) # makes a few transformations easier akin to dplyr and tidyr
library(mgcv) # needed for fitting GAMs
library(broom) # helps with tidying output of statistical models for further processing
library(ggplot2) # for plots
library(grid) # for copying text from console
library(pheatmap)

#pull data from KEGG to use for gene sets
hsa_kegg = download_KEGG("hsa")

# so the KEGG pathway ID to entrez ID is what we want but gsva wants it in list format 
# so merge the columns and split them so each pathway in the list is a a vector containing the entrez IDs
kegg_map = merge(hsa_kegg$KEGGPATHID2EXTID, hsa_kegg$KEGGPATHID2NAME,
                  by.x = "from", by.y = "from")
gsets = split(as.character(kegg_map$to.x), kegg_map$to.y)

# Variance Stabilising Transformation (VST) for GSVA:
vsd = vst(ddsTC, blind = FALSE)
head(assay(vsd), 3)
colData(vsd)

# Take the expression set
expr_vst = assay(vsd)
# check all samples are in the metadata
stopifnot(all(colnames(expr_vst) == meta$ID))

# pull the uniprot IDs
prot_names = read.table("metadata/gene2uniprotID.tsv", header=F, sep="\t")
#rownames(prot_names) = prot_names[,1]
prots = prot_names[ -c(3:4) ]
colnames(prots) = c("geneID", "uniprot_gn_id")

# There are 18499 genes identified in counts (some of which will be non-coding), 14720 genes identified with a uniprot ID
# 8594 of the uniprot IDs have mapped to a unique human/non-human homolog (meaning there are 6126 duplicates)
# Some of these could be poor mapping and some will be paralogs or gene variants which map to the same human protein.
nrow(vsd)
length(prots$uniprot_gn_id)
length(unique(prots$uniprot_gn_id))

# The duplicates will not be appropriate for some gene set analyses and cause issue with scoring. 
# GSEA assumes gene independence for example. GSVA doesn't but I don't have gene sets and need to obtain them from somewhere, also pathview and KEGG require unique gene IDs.
# I will choose to keep the gene with the highest expression mean to try and keep the best representation
gene_means = rowMeans(expr_vst)  
gene_means = data.frame(geneID = names(gene_means), meanExpr = gene_means)
prots = merge(prots, gene_means, by="geneID")

prots_filtered = prots %>%
  group_by(uniprot_gn_id) %>%
  slice_max(order_by = meanExpr, n = 1) %>%
  ungroup()

# Duplicates sucessfully removed :>
length(prots_filtered$uniprot_gn_id)
length(unique(prots_filtered$uniprot_gn_id))

# Convert the uniprot IDs to entrez IDs for use cases where they are required (such as pathview because of KEGG)
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
map = getBM(
  attributes = c("uniprot_gn_id", "entrezgene_id"),
  filters    = "uniprot_gn_id",
  values     = prots_filtered$uniprot_gn_id,
  mart       = mart
)


# of the 14720 proteins identified with a uniprot ID, and after filtering for unique mappings we then see 7543 genes
# 7542 of these are unique (so only 1 left is a duplicate at this point, most of the duplication was taken out with the duplicate uniprots)
length(map$entrezgene_id)
length(unique(map$entrezgene_id))

prots_filtered = merge(prots_filtered, map, by = "uniprot_gn_id")

prots_twice_filt = prots_filtered %>%
  group_by(entrezgene_id) %>%
  slice_max(order_by = meanExpr, n = 1) %>%
  ungroup()

# and we finally have ONLY unique gene IDs
length(prots_twice_filt$entrezgene_id)
length(unique(prots_twice_filt$entrezgene_id))

# there are multiple entrezIDs per gene ID though so I also need to remove those, I don't have a base mean to go on here so I'll just remove the first each time
length(prots_twice_filt$geneID)
length(unique(prots_twice_filt$geneID))

prots_final = prots_twice_filt %>%
  group_by(uniprot_gn_id) %>%
  slice(1) %>%
  ungroup()

# finally we have 7475 unique entrez IDs to apply to our expr_vst
length(prots_final$entrezgene_id)
length(unique(prots_final$entrezgene_id))
length(prots_final$geneID)
length(unique(prots_final$geneID))

prots_final = data.frame(prots_final)
prots_final = prots_final[ -c(1,3) ]
prots_final = na.omit(prots_final)
rownames(prots_final) = prots_final[,1]
expr_vst = merge(expr_vst, prots_final, by = "row.names")
expr_vst = expr_vst %>% dplyr::select(-geneID)
row.names(expr_vst)= expr_vst$entrezgene_id
expr_vst = expr_vst %>% dplyr::select(-Row.names, -entrezgene_id)

# save prots_final for later kegg use
write.table(prots_final, "sets/kegg_gene2entrez.tsv", row.names= F, sep ="\t", quote = F)

# Finally, make sure the expr_vst is a matrix again
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

# File name and the indicators are messing up the pivotting of the dataframe of results 
meta_clean = meta %>% dplyr::select(-c(FileName, Azo, Cyp, Imid, Pro))

# make a long format dataframe of results
df = as.data.frame(t(gsva_res)) |>
  rownames_to_column("ID") |>
  left_join(meta_clean, by="ID") |>
  pivot_longer(cols = -c(ID,Treatment,Time,Timetreatment,Rep),
               names_to="Pathway", values_to="Score")

df$Time = as.numeric(as.character(df$Time))
df$Rep = as.factor(df$Rep)
df$Treatment = relevel(df$Treatment, ref = "Control")

fit_gam_one = function(d){
  # model GAMs with treatment-specific smoothing + main treatment effect over time 
  # setting 6 as maximum knots (6 time points - which can be computationally intensive but won't sacrifice accuracy) 
  # and using  maximum likelihoods
  gam(Score ~ Treatment + s(Time, by=Treatment, k=6), data=d, method="ML")
}
# fit a model per pathway (so each treatment will have a spline per KEGG pathway - 
# later we will both visualise this and compare to a treatment averaged spline per treatment and control per pathway)
# essentially doing a LRT like in the overall deseq2 model but this time on each treatment with the same model 
# but then making each treatment share the spline with the control, and seeing if the spline trajectory is significantly different
paths = unique(df$Pathway)
model_fits  = setNames(vector("list", length(paths)), paths)
for(p in paths){
  model_fits[[p]] = tryCatch(fit_gam_one(filter(df, Pathway==p)))
}

# Need to check the model fits to ensure they are actually good, representative hits and not just cramming data in a model and throwing it at the wall
dir.create("GAM_diagnostics")
dir.create("GAM_diagnostics/Treatments/")
dir.create("GAM_diagnostics/Treatments/KEGG")

for(pathway in names(model_fits)) {
  fit = model_fits[[pathway]]
  
  if(inherits(fit, "try-error") || is.null(fit)) next  # skip failed models
  
  safe_pathway = gsub("[^A-Za-z0-9]", "_", pathway)
  
  pdf(paste0("GAM_diagnostics/Treatments/KEGG/", safe_pathway ,".pdf"), width=10, height=7)
  
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
flags = data.frame(Pathway=names(model_fits), PoorFit=NA)
for(i in seq_along(model_fits)) {
  fit = model_fits[[i]]
  if(inherits(fit, "try-error") || is.null(fit)) next
  k_table = mgcv:::k.check(fit)
  res = residuals(fit, type="deviance")
  
  flags$PoorFit[i] = any(k_table[,"p-value"] < 0.05, abs(res) > 3)
}
bad_fits = subset(flags, flags$PoorFit == T)

# There are 9 pathways modelled to have a potentially bad fit, so I will inspect their plots individually :>
# Coronavirus disease - COVID-19 not looking too good, showing potential signs of heteroskedasticity. Will drop this...
# Ribosome also showing signs of heteroskedasticity
# Virion - Ebolavirus, Lyssavirus and Morbillivirus also shows signs of heteroskedasticity. 
# Rest are fine - will not consider these 3 as valid results.

#write.csv(bad_fits, "GAM_diagnostics/poor_fits.csv", row.names=FALSE)

# compute time-change significance from fitted GAMs (aka if the GSVA score changes signifcantly)
time_sig_res = lapply(names(model_fits), function(pathway) {
  fit = model_fits[[pathway]]
  if(inherits(fit, "try-error") || is.null(fit)) return(NULL)
  
  s_summ = summary(fit)$s.table  # edf, Ref.df, F, p-value
  
  tibble(
    Pathway = pathway,
    Treatment = gsub("s\\(Time\\):Treatment", "", rownames(s_summ)),  # remove the prefix
    p_value_time = s_summ[, "p-value"],
    time_sig_flag = ifelse(s_summ[, "p-value"] < 0.05, "T", "F")
  )
})
time_sig_res = bind_rows(time_sig_res)


#########
# function to compare each treatment against control and the treatment for a given pathway (nested model to perform LRT on)
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

# run all f the control comparison tests pairwise across all treatments with the control within each pathway (this takes a good min)
pairwise_res = lapply(paths, function(p){
  d = filter(df, Pathway==p)
  tmp = compare_to_control(d, control="Control") |>
    mutate(Pathway=p)
})

# bind the rows, do FDR mitigation and 
# make a flag for if the results for the treatment that were significantly improved from being modelled separately to the control or not which will be helpful for plotting
pairwise_res = bind_rows(pairwise_res) |>
  mutate(p_adj = p.adjust(p_value, "BH"),
         sig_flag = ifelse(p_adj < 0.05, "T", "F"))

# add the significance flags into the gsva dataframe
df_with_sig = df %>%
  left_join(pairwise_res %>% dplyr::select(Pathway, Treatment, sig_flag),
            by = c("Pathway", "Treatment"))

## aslo add the sig flags for trajectory changing over time to the plot :>
df_with_sig = df_with_sig %>%
  left_join(time_sig_res %>% dplyr::select(Pathway, Treatment, time_sig_flag),
            by = c("Pathway", "Treatment"))

df_with_sig = df_with_sig %>%
  filter(!Pathway %in% c("Coronavirus disease - COVID-19", "Ribosome", "Virion - Ebolavirus, Lyssavirus and Morbillivirus"))

################# PLOTTING ##############################

df_with_sig = df_with_sig %>%
  mutate(Treatment = factor(Treatment, levels = c("Control", setdiff(unique(Treatment), "Control"))))

# Define treatment colors (adjust as needed)
treatment_colors = c(
  "Control" = "black", "Azo-Pro" = "#AE4D66", "Azoxy" = "#FF9933",
  "Cyp" = "#00CC33", "Cyp-Pro" = "#2E6666", "Imd-Pro" = "#4866CC",
  "Imid" = "#33CCFF", "Pro" = "#5C0099"
)

treatment_sets = list(
  set1 = c("Control", "Pro", "Azoxy", "Azo-Pro"),
  set2 = c("Control", "Pro", "Cyp", "Cyp-Pro"),
  set3 = c("Control", "Pro", "Imid", "Imd-Pro")
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
dir.create("figures/GSVA_Treatment_Splines/KEGG_pathways")

# Define treatment sets
treatment_sets = list(
  "Azo-Pro Set" = c("Control", "Pro", "Azoxy", "Azo-Pro"),
  "Cyp-Pro Set" = c("Control", "Pro", "Cyp", "Cyp-Pro"),
  "Imd-Pro Set" = c("Control", "Pro", "Imid", "Imd-Pro")
)

# Loop over pathways

for (pathway in unique(df_with_sig$Pathway)) {
  
  # Subset data for this pathway
  d = df_with_sig %>% filter(Pathway == pathway)
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
  ggsave(paste0("figures/GSVA_Treatment_Splines/KEGG_pathways/", safe_pathway,".png"), p,
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

# fit a model per Pathway so each predictor and interaction term will have a spline per KEGG_Pathways - 
# later we will visualise this, although mostly just the interaction splines (above zero shows synergy and below zero shows antagonism)
# will also do essentially an LRT with a reduced model with no interaction terms again, to see if it shows significant model improvement (aka presence of non additive interactions)
paths = unique(df$Pathway)
int_model_fits  = setNames(vector("list", length(paths)), paths)
for(p in paths){
  int_model_fits[[p]] = tryCatch(fit_int_gam(filter(df, Pathway==p)))
}

# Need to check the model fits to ensure they are actually good, representative hits and not just cramming data in a model and throwing it at the wall
dir.create("GAM_diagnostics")
dir.create("GAM_diagnostics/Interactions/")
dir.create("GAM_diagnostics/Interactions/KEGG")

for(pathway in names(int_model_fits)) {
  int_fit = int_model_fits[[pathway]]
  
  if(inherits(int_fit, "try-error") || is.null(int_fit)) next  # skip failed models
  
  safe_pathway = gsub("[^A-Za-z0-9]", "_", pathway)
  
  pdf(paste0("GAM_diagnostics/Interactions/KEGG/", safe_pathway ,".pdf"), width=10, height=7)
  
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
int_flags = data.frame(Pathway=names(int_model_fits), PoorFit=NA)
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
  int_d = filter(df, Pathway==p)
  int_tmp = reduced_interactions(int_d) |>
    mutate(Pathway=p)
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
  left_join(int_pairwise_res %>% dplyr::select(Pathway, Interaction, sig_flag),
            by = c("Pathway", "Interaction"))


#### Plotting Interactions :> ####

dir.create("figures/GSVA_Interaction_Splines")
dir.create("figures/GSVA_Interaction_Splines/KEGG_Pathways")

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

for (pathway in unique(int_df_with_sig$Pathway)) {
  # Subset data for this pathway
  d = int_df_with_sig %>% filter(Pathway == pathway)
  fit = int_model_fits[[pathway]]
  if (inherits(fit, "try-error") || is.null(fit)) next  # skip failed models
  
  safe_pathway = gsub("[^A-Za-z0-9]", "_", pathway)
  
  # Grab p-values for this pathway
  sig_res = int_pairwise_res %>% bind_rows() %>%
    filter(Pathway == pathway) %>%
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
  
  ggsave(paste0("figures/GSVA_Interaction_Splines/KEGG_Pathways/", safe_pathway,".png"), p,
         width = 9, height = 6)
}


#### List genes from KEGG sets that are also in prot_final for later working out overlap with GO ####

# prots_final contains the Entrez IDs of the genes actually measured
possible_genes = prots_final$entrezgene_id

# filtering each KEGG pathway to keep only genes in universe
kegg_sets_filtered = lapply(gsets, function(gs) {
  intersect(gs, possible_genes)
})

saveRDS(kegg_sets_filtered, file = "sets/kegg_sets_entrez.rds")

# make a vector to rename the entrez IDs back for comparison to go.
entrez_to_geneID = setNames(prots_final$geneID, prots_final$entrezgene_id)

# map them back
kegg_sets_final = lapply(kegg_sets_filtered, function(genes_entrez) {
  # Replace Entrez IDs with geneIDs
  geneIDs = entrez_to_geneID[genes_entrez]
  # Remove any NAs (in case some Entrez IDs are not found)
  geneIDs = geneIDs[!is.na(geneIDs)]
  return(geneIDs)
})

dir.create("sets")

saveRDS(kegg_sets_final, file = "sets/kegg_sets_final.rds")

#### save project ####

# Save .RData file
save.image(file=paste0("KEGG_GSVA_Results", ".RData"))
