## Script to filter DESeq2 output and create gene results tables both from the LRT and from each individual Wald test.

# This DOES NOT create DEG lists, only filters down using the omnibus test method and produces results tables with all relevant FDR/FWER controlled p values (AND includes the L2FC for later)
# Again, this script takes a few hours to run with my data so it's advisable to submit the script as a slurm job or run in the background with a nohup command - don't try to run it actively and wait.

#### Set-up ####
rm(list=ls())
#Set working directory
workDir = "/mnt/clusters/admiral/data/c21082179/BIT104/New_Analysis"  # working directory for the R session
setwd(workDir)


## Read in the .RData file
load("DESeq2_TC_Results.RData")

##libraries
library(pheatmap)
library(DESeq2)
library("stageR")



#### IMPORTANT TO NOTE ####
# I will initially be taking p VALUES NOT p-ADJUSTED values from the LRT1 and soon to be performed Wald tests, as these are going to be corrected for AFTERWARD in a 2-stage FDR correction with stageR
# Otherwise ALWAYS take p-adjusted vals to mitigate inflated false positives.



#### 1. TIME-AVERAGED INDIVIDUAL CHEMICAL ANALYSIS ####

# For the treatment only I'm going to essentially take an average, which will require messing about with the contrast argument see: https://github.com/tavareshugo/tutorial_DESeq2_contrasts/blob/main/DESeq2_contrasts.md
# that's the concept I'm using - he has since adapted that to warn against using it for partially crossed or imbalanced designed but mine is fully crossed
# I do have more control replicates but that is not directly affecting anything here as ALL timepoints within a treatment have the same amount of replicates
# so I'm not going to mess about trying to do weighting - but for future reference weighting needs to be manually adapted if there are different replicate amounts.

# Assign the UNadjusted p values from LRT1 to be the screening stage 1 for the stageR FDR correction, and provide the new object with the gene names
p_screen = resLRT1$pvalue
names(p_screen) = rownames(resLRT1)

# Get the model matrix to use for the contrast vector
mod_mat = model.matrix(design(WLD), colData(WLD))

# pull different treatment types (Except control)
treatments = setdiff(levels(WLD$Treatment), "Control")

#lists for loop
p_list = list()
log2FC_list = list()
log2FC_shrk_list = list()

for(tr in treatments) {
  
  # here I am basically pullihng all the coefficients associated with the treatment (aka the treatment at each timepoint)
  # the main effect is just the treatment at timepoint 1h, basically this is considered the baseline effect and all the other effects add onto it :>
  main_col = paste0("Treatment", tr)
  
  # get all the other effects/coefficients (aka everything at the other timepoints, which are each the baseline plus additional expression at that timepoint)
  inter_cols = grep(paste0("Treatment", tr, ":Time"), colnames(mod_mat), value = TRUE)
  
  # make an empty contrast vector, it needs the names of the model matrix and everything should be 0 at default before adding the weights
  contrast_vec = rep(0, ncol(mod_mat))
  names(contrast_vec) = colnames(mod_mat)
  
  # combine the coefficient columns so we have the baseline and the interaction effects, then weight them because the're supposed to add to 1
  all_cols = c(main_col, inter_cols)
  contrast_vec[all_cols] = 1 / length(all_cols)
  
  # extract results using numeric contrast vector, and get the shrunken L2FCs because they're more reliable/conservative
  res = results(WLD, contrast = contrast_vec, alpha = 0.05)
  res_shrk = lfcShrink(WLD, contrast = contrast_vec, type = "ashr", res = res)
  
  # put in lists because stageR needs stuff in lists
  coef_name = paste0(tr, "_vs_control")
  p_list[[coef_name]] = res$pvalue
  log2FC_list[[paste0(coef_name, "_L2FC")]] = res$log2FoldChange
  log2FC_shrk_list[[paste0(coef_name, "_shrk_L2FC")]] = res_shrk$log2FoldChange
}

# Restructuring to prepare for stage R
p_confirm = do.call(cbind, p_list)
rownames(p_confirm) = rownames(WLD)
l2fc_confirm = do.call(cbind, log2FC_list)
rownames(l2fc_confirm) = rownames(WLD)
l2fc_shrk_confirm = do.call(cbind, log2FC_shrk_list)
rownames(l2fc_shrk_confirm) = rownames(WLD)


# Now a stage R object will be created which will contain the information required for stageR analysis to be run in the format required (including the omnibus test results: pscreen, and the wald test results: pconfirm)
# I have also clarified to th stage R function that the p values for the LRT were NOT adjusted prior
stg_r_ob = stageR(
  pScreen = p_screen,
  pConfirmation = p_confirm,
  pScreenAdjusted = F
)

# Here we are performing the stagewise adjustment, this tests ACROSS all genes (FDR corrected with BH) for significance p.adj <0.05
# Genes which PASS this screening stage show significant DE associated with the treatment variable.
# Any genes that pass the screening stage are tested at the confirmation stage where FWER or FDR are controlled ACROSS hypotheses rather than genes this time
# This looks at each treatment comparison individually, and in this case we control the within-gene family wise error rate with Holm, alpha 0.05.
# This is a very stringent, conservative control of FWER, and controls the probability of making at least 1 false positive across all hypotheses.
# I'm going to do a version of this with BH later which will be used for pathways and exploratory analysis where we want to be less conservative, but the Holm method is great for when we want to look at specific genes/tests rather than pathways.
stg_r_ob = stageWiseAdjustment(stg_r_ob, method = "holm", alpha = 0.05)


# extract results
stage_r_res = getAdjustedPValues (stg_r_ob, order = F, onlySignificantGenes = F)

# Note for my output log
paste("Check the number of DE instances between treatments identified (Holm method):")
paste("WARNING: These are not necessarily statistically detectable in this analysis - that will depend on the LFC, which will be detected in the DEG_list_maker.R script")
sum(stage_r_res[ , -1] < 0.05, na.rm = TRUE)

# Note for output log
paste("Check number of genes where at least one DE instance between treatment was identified (Holm method):")
paste("WARNING: These are not necessarily statistically detectable in this analysis - that will depend on the LFC, which will be detected in the DEG_list_maker.R script")
sum(apply(stage_r_res[ , -1], 1, function(row) any(row < 0.05, na.rm = TRUE)))

## Note for output log - this is hashed out because I don't actually have these to compare now I have combined controls.
#paste("Check instances where carrier and water show DE (Holm method):")
#sum(stage_r_res[ , "Carrier_vs_Water"] < 0.05, na.rm = TRUE)

# Add the log 2 fold change values for the wald tests to the results table
res1_tmp = merge(stage_r_res, l2fc_confirm, by = "row.names")
rownames(res1_tmp)=res1_tmp$Row.names
res1_tmp$Row.names = NULL
res1_table = merge(res1_tmp, l2fc_shrk_confirm, by = "row.names")
rownames(res1_table)=res1_table$Row.names
res1_table$Row.names = NULL

# Save results to a table
write.table(res1_table, "DEG_lists/res1_table.tsv", row.names= T, sep ="\t", quote = F)


## This is now going to be exactly the same thing but I will use BH - this isn't explicitly provided as a method WITHIN stageR, so I have to apply it to pconfirm outside of stage R and then select method="none".
# BH is less conservative it controls the PROPORTION of false positives and thus will be better for identifying DE pathways which could be more negatively affected by false negatives
BH_adjusted_matrix = t(apply(p_confirm, 1, function(p) p.adjust(p, method = "BH")))

stage_r_BH = stageR(pScreen = p_screen,
                    pConfirmation = BH_adjusted_matrix,
                    pScreenAdjusted = FALSE)

stage_r_BH = stageWiseAdjustment(object = stage_r_BH,
                                   method = "none", # Using none here because I have already BH adjusted in this case
                                   alpha = 0.05)

# Extract results
stage_r_BH_res = getAdjustedPValues (stage_r_BH, order = F, onlySignificantGenes = F)

# Note for output log
paste("Check the number of DE instances between treatments identified (BH method):")
paste("WARNING: These are not necessarily statistically detectable in this analysis - that will depend on the LFC, which will be detected in the DEG_list_maker.R script")
sum(stage_r_BH_res[ , -1] < 0.05, na.rm = TRUE)

# Note for output log
paste("Check number of genes where at least one DE instance between treatment was identified (BH method):")
paste("WARNING: These are not necessarily statistically detectable in this analysis - that will depend on the LFC, which will be detected in the DEG_list_maker.R script")
sum(apply(stage_r_BH_res[ , -1], 1, function(row) any(row < 0.05, na.rm = TRUE)))

## Note for output log (again don't need this now I have combined controls but this was important for me to look at prior to determine if I COULD combine them)
#paste("Check instances where carrier and water show DE (BH method):")
#sum(stage_r_BH_res[ , "Carrier_vs_Water"] < 0.05, na.rm = TRUE)

# Add the log 2 fold change values for the wald tests to the results table
res1_BH_tmp = merge(stage_r_BH_res, l2fc_confirm, by = "row.names")
rownames(res1_BH_tmp)=res1_BH_tmp$Row.names
res1_BH_tmp$Row.names = NULL
res1_BH_table = merge(res1_BH_tmp, l2fc_shrk_confirm, by = "row.names")
rownames(res1_BH_table)=res1_BH_table$Row.names
res1_BH_table$Row.names = NULL

# Save results
write.table(res1_BH_table, "DEG_lists/res1_BH_table.tsv", row.names= T, sep ="\t", quote = F)






#### 2. INDIVIDUAL CHEMICAL PER TIME POINT ANALYSIS ####
# This version is for GSEA analysis over time, aka the middle ground: more detailed than the treatment overview but focuses on all treatment associated DE at each timepoint
# so uses LRT1 again to filter and not quite as zoomed in and sensitive to treatment jumps over time as the within treatment time comparisons are...

# Assign the UNadjusted p values from LRT1 to be the screening stage 1 for the stageR FDR correction, and provide the new object with the gene names
p_screen = resLRT1$pvalue
names(p_screen) = rownames(resLRT1)

# pull different treatment types (Except control) and timepoint types
treatments = setdiff(levels(WLD$Treatment), "Control")
timepoints = levels(WLD$Time)

treatments_dots = gsub("-", ".", treatments)


#lists for loop
p_list = list()
log2FC_list = list()
log2FC_shrk_list = list()


for (i in 1:length(treatments)) {
  treatment = treatments[i]
  treatment_dots = treatments_dots[i]
  
  for (time in timepoints) {
    
    contrast_label = paste(treatment, time, "vs_Control", time, sep = "_")
    
    # build coefficient names exactly as in resultsNames(WLD)
    coef1 = paste0("Time_", time, "_vs_1")
    coef2 = paste0("Treatment", treatment_dots, ".Time", time)
    coef_baseline = paste0("Treatment_", treatment_dots, "_vs_Control")
    
    if (time == "1") {
      
      # When time is baseline, only use Treatment
      if (coef_baseline %in% resultsNames(WLD)) {
        res = results(WLD, name = coef_baseline, alpha = 0.05)
        res_shrk = lfcShrink(WLD, coef = coef_baseline, type="ashr", res=res)
        
        p_list[[contrast_label]] = res$pvalue
        log2FC_list[[paste0(contrast_label, "_L2FC")]] = res$log2FoldChange
        log2FC_shrk_list[[paste0(contrast_label, "_shrk_L2FC")]] = res_shrk$log2FoldChange
        
      } else {
        message(paste("Skipping", contrast_label, "- coefficient not found"))
      }
    } else {
      # time point are not baseline: contrast of two coefficients
      
      if (all(c(coef1, coef2) %in% resultsNames(WLD))) {
        res = results(WLD, contrast = list(c(coef_baseline, coef1, coef2)), alpha = 0.05)
        res_shrk = lfcShrink(WLD, contrast = list(c(coef_baseline, coef1, coef2)), type="ashr", res=res)
        
        p_list[[contrast_label]] = res$pvalue
        log2FC_list[[paste0(contrast_label, "_L2FC")]] = res$log2FoldChange
        log2FC_shrk_list[[paste0(contrast_label, "_shrk_L2FC")]] = res_shrk$log2FoldChange
      } else {
        message(paste("Skipping", contrast_label, "- coefficients not found"))
      }
    }
  }
}




p_confirm = do.call(cbind, p_list)
rownames(p_confirm) = rownames(WLD)
l2fc_confirm = do.call(cbind, log2FC_list)
rownames(l2fc_confirm) = rownames(WLD)
l2fc_shrk_confirm = do.call(cbind, log2FC_shrk_list)
rownames(l2fc_shrk_confirm) = rownames(WLD)

stg_r_ob = stageR(
  pScreen = p_screen,
  pConfirmation = p_confirm,
  pScreenAdjusted = F
)

stg_r_ob = stageWiseAdjustment(stg_r_ob, method = "holm", alpha = 0.05)


# Extract results
stage_r_res = getAdjustedPValues (stg_r_ob, order = F, onlySignificantGenes = F)

# Note for output log
paste("Check the number of DE instances between treatments and control at a certain timepoint (Holm method):")
paste("WARNING: These are not necessarily statistically detectable in this analysis - that will depend on the LFC, which will be detected in the DEG_list_maker.R script")
sum(stage_r_res[ , -1] < 0.05, na.rm = TRUE)

# Note for output log
paste("Check number of genes where at least one DE instance between treatments and control at a certain timepoint was identified (Holm method):")
paste("WARNING: These are not necessarily statistically detectable in this analysis - that will depend on the LFC, which will be detected in the DEG_list_maker.R script")
sum(apply(stage_r_res[ , -1], 1, function(row) any(row < 0.05, na.rm = TRUE)))

# Add the log 2 fold change values for the wald tests to the results table
res2_tmp = merge(stage_r_res, l2fc_confirm, by = "row.names")
rownames(res2_tmp)=res2_tmp$Row.names
res2_tmp$Row.names = NULL
res2_table = merge(res2_tmp, l2fc_shrk_confirm, by = "row.names")
rownames(res2_table)=res2_table$Row.names
res2_table$Row.names = NULL

# Save results
write.table(res2_table, "DEG_lists/res2_table.tsv", row.names= T, sep ="\t", quote = F)


## BH edit
BH_adjusted_matrix = t(apply(p_confirm, 1, function(p) p.adjust(p, method = "BH")))

stage_r_BH = stageR(pScreen = p_screen,
                    pConfirmation = BH_adjusted_matrix,
                    pScreenAdjusted = FALSE)

stage_r_BH = stageWiseAdjustment(object = stage_r_BH,
                                 method = "none", # Using none here because I have already BH adjusted in this case
                                 alpha = 0.05)

# Extract results
stage_r_BH_res = getAdjustedPValues (stage_r_BH, order = F, onlySignificantGenes = F)

# Note for output log
paste("Check the number of DE instances between treatments and control at different timepoints (BH method):")
paste("WARNING: These are not necessarily statistically detectable in this analysis - that will depend on the LFC, which will be detected in the DEG_list_maker.R script")
sum(stage_r_BH_res[ , -1] < 0.05, na.rm = TRUE)

# Note for output log
paste("Check number of genes where at least one DE instance between treatments and control at a certain timepoint was identified (BH method):")
paste("WARNING: These are not necessarily statistically detectable in this analysis - that will depend on the LFC, which will be detected in the DEG_list_maker.R script")
sum(apply(stage_r_BH_res[ , -1], 1, function(row) any(row < 0.05, na.rm = TRUE)))

# Add the log 2 fold change values for the wald tests to the results table
res2_BH_tmp = merge(stage_r_BH_res, l2fc_confirm, by = "row.names")
rownames(res2_BH_tmp)=res2_BH_tmp$Row.names
res2_BH_tmp$Row.names = NULL
res2_BH_table = merge(res2_BH_tmp, l2fc_shrk_confirm, by = "row.names")
rownames(res2_BH_table)=res2_BH_table$Row.names
res2_BH_table$Row.names = NULL

# Save results
write.table(res2_BH_table, "DEG_lists/res2_BH_table.tsv", row.names= T, sep ="\t", quote = F)




#### 3. TIME DYNAMICS WITHIN CHEMICALS ANALYSIS ####

# Here what I'm doing is very similar to the previous method but here genes which pass this screening stage show significant DE associated with the treatment:time interaction term.
# This means that they show both differential expression due to time-dependent treatment affects, which will identify interesting treatment based behaviour of gene expression over time.
# We will then be able to compare individual time points WITHIN each treatment.
# This contrasting method is messy and complicated :< but I am fairly confident I actually have it correct now
p_screen = resLRT2$pvalue
names(p_screen) = rownames(resLRT2)

# Loop through all of the TREATMENT Wald tests to be performed and extract the p values and log fold changes for stage 2 FWER correction.
treatment_coefs = levels(WLD$Treatment)[-1] # This time excluding the treatment control because we're not looking at things BETWEEN the treatments, we're looking at differences between time points WITHIN the treatments as opposed to between time points in the control.
## TO CLARIFY: this is not excluding the control from calculations... the control was already used as a BASELINE to determine the results of the LRT2 showing association with treatment:time
# essentially we are looking at changes between time points in the treatments AS COMPARED to time points in the control - so the control won't show anything
#We can always visually compare interesting finds to the control later, we just don't need to investigate the control itself here.

# They're also formatted a bit dodgily in the treatment time results names so I need to edit that and replace all hyphens with periods.
treatment_dots_coefs = gsub("-", ".", treatment_coefs)

timepoints = levels(WLD$Time)
time_pairs = combn(timepoints, 2, simplify = FALSE)

p_list = list()
log2FC_list = list()
log2FC_shrk_list = list()

for (i in 1:length(treatment_coefs)) {
  treatment = treatment_coefs[i]
  treatment_dots = treatment_dots_coefs[i]
  for (pair in time_pairs) {
    time1 = pair[1]
    time2 = pair[2]
    contrast_label = paste(treatment, "time", time2, "vs", time1, sep = "_")
    
    # interaction coefficients
    coef1 = paste0("Treatment", treatment_dots, ".Time", time2)  # interaction at time2
    coef2 = paste0("Treatment", treatment_dots, ".Time", time1)  # interaction at time1
    
    # time main-effect coefficients (control time effects)
    time1_coef = paste0("Time_", time1, "_vs_1")
    time2_coef = paste0("Time_", time2, "_vs_1")
    
    if (time1 == timepoints[1]) {   # baseline = first level, e.g. "1"
      # DrugA_time2 vs DrugA_time1(=baseline) = Time_t2 + TreatmentA.Time_t2
      if (all(c(time2_coef, coef1) %in% resultsNames(WLD))) {
        res = results(WLD, contrast = list(c(time2_coef, coef1)), alpha = 0.05)
        res_shrk = lfcShrink(WLD, contrast = list(c(time2_coef, coef1)), type="ashr", res=res)
        p_list[[contrast_label]] = res$pvalue
        log2FC_list[[paste0(contrast_label, "_L2FC")]] = res$log2FoldChange
        log2FC_shrk_list[[paste0(contrast_label, "_shrk_L2FC")]] = res_shrk$log2FoldChange
      } else {
        message("Skipping ", contrast_label, " - needed coefficients missing")
      }
    } else {
      # General case: (Time_t2 + T:Time_t2) - (Time_t1 + T:Time_t1)
      if (all(c(time2_coef, coef1, time1_coef, coef2) %in% resultsNames(WLD))) {
        res = results(WLD, contrast = list(c(time2_coef, coef1), c(time1_coef, coef2)), alpha = 0.05)
        res_shrk = lfcShrink(WLD, contrast = list(c(time2_coef, coef1), c(time1_coef, coef2)), type="ashr", res=res)
        p_list[[contrast_label]] = res$pvalue
        log2FC_list[[paste0(contrast_label, "_L2FC")]] = res$log2FoldChange
        log2FC_shrk_list[[paste0(contrast_label, "_shrk_L2FC")]] = res_shrk$log2FoldChange
      } else {
        message("Skipping ", contrast_label, " - needed coefficients missing")
      }
    }
  }
}


p_confirm = do.call(cbind, p_list)
rownames(p_confirm) = rownames(WLD)
l2fc_confirm = do.call(cbind, log2FC_list)
rownames(l2fc_confirm) = rownames(WLD)
l2fc_shrk_confirm = do.call(cbind, log2FC_shrk_list)
rownames(l2fc_shrk_confirm) = rownames(WLD)

stg_r_ob = stageR(
  pScreen = p_screen,
  pConfirmation = p_confirm,
  pScreenAdjusted = F
)

stg_r_ob = stageWiseAdjustment(stg_r_ob, method = "holm", alpha = 0.05)


# Extract results
stage_r_res = getAdjustedPValues (stg_r_ob, order = F, onlySignificantGenes = F)

# Note for output log
paste("Check the number of DE instances between time points in different treatments compared to the same timepoints in the control identified (Holm method):")
paste("WARNING: These are not necessarily statistically detectable in this analysis - that will depend on the LFC, which will be detected in the DEG_list_maker.R script")
sum(stage_r_res[ , -1] < 0.05, na.rm = TRUE)

# Note for output log
paste("Check number of genes where at least one DE instance between time points in a treatment was identified (Holm method):")
paste("WARNING: These are not necessarily statistically detectable in this analysis - that will depend on the LFC, which will be detected in the DEG_list_maker.R script")
sum(apply(stage_r_res[ , -1], 1, function(row) any(row < 0.05, na.rm = TRUE)))

# Add the log 2 fold change values for the wald tests to the results table
res3_tmp = merge(stage_r_res, l2fc_confirm, by = "row.names")
rownames(res3_tmp)=res3_tmp$Row.names
res3_tmp$Row.names = NULL
res3_table = merge(res3_tmp, l2fc_shrk_confirm, by = "row.names")
rownames(res3_table)=res3_table$Row.names
res3_table$Row.names = NULL

# Save results
write.table(res3_table, "DEG_lists/res3_table.tsv", row.names= T, sep ="\t", quote = F)


## BH edit
BH_adjusted_matrix = t(apply(p_confirm, 1, function(p) p.adjust(p, method = "BH")))

stage_r_BH = stageR(pScreen = p_screen,
                     pConfirmation = BH_adjusted_matrix,
                     pScreenAdjusted = FALSE)

stage_r_BH = stageWiseAdjustment(object = stage_r_BH,
                                  method = "none", # Using none here because I have already BH adjusted in this case
                                  alpha = 0.05)

# Extract results
stage_r_BH_res = getAdjustedPValues (stage_r_BH, order = F, onlySignificantGenes = F)

# Note for output log
paste("Check the number of DE instances between times within treatments identified (BH method):")
paste("WARNING: These are not necessarily statistically detectable in this analysis - that will depend on the LFC, which will be detected in the DEG_list_maker.R script")
sum(stage_r_BH_res[ , -1] < 0.05, na.rm = TRUE)

# Note for output log
paste("Check number of genes where at least one DE instance between times within treatments was identified (BH method):")
paste("WARNING: These are not necessarily statistically detectable in this analysis - that will depend on the LFC, which will be detected in the DEG_list_maker.R script")
sum(apply(stage_r_BH_res[ , -1], 1, function(row) any(row < 0.05, na.rm = TRUE)))

# Add the log 2 fold change values for the wald tests to the results table
res3_BH_tmp = merge(stage_r_BH_res, l2fc_confirm, by = "row.names")
rownames(res3_BH_tmp)=res3_BH_tmp$Row.names
res3_BH_tmp$Row.names = NULL
res3_BH_table = merge(res3_BH_tmp, l2fc_shrk_confirm, by = "row.names")
rownames(res3_BH_table)=res3_BH_table$Row.names
res3_BH_table$Row.names = NULL

# Save results
write.table(res3_BH_table, "DEG_lists/res3_BH_table.tsv", row.names= T, sep ="\t", quote = F)






#### 4. TIME-AVERAGED CHEMICAL INTERACTION ANALYSIS ####

# Assign the UNadjusted p values from LRT1 to be the screening stage 1 for the stageR FDR correction, and provide the new object with the gene names
p_screen = resINTLRT1$pvalue
names(p_screen) = rownames(resINTLRT1)

# Get the model matrix to use for the contrast vector
mod_mat = model.matrix(design(INTWLD), colData(INTWLD))

# identify different interaction types
interactions = c("Azo_Pro", "Cyp_Pro", "Imid_Pro")

#lists for loop
p_list = list()
log2FC_list = list()
log2FC_shrk_list = list()

for(int in interactions) {
  
  # Extract the changing chem
  i = strsplit(int, "_Pro")[[1]]
  
  # here I am basically pulling all the coefficients associated with the treatment (aka the treatment at each timepoint) as seen in the mod_mat
  first_col = paste0(i, "pres:Propres")
  
  # get all the other effects/coefficients (aka everything at the other timepoints, which are each the baseline plus additional expression at that timepoint)
  other_cols = grep(paste0(first_col, ":Time"), colnames(mod_mat), value = TRUE)
  
  # make an empty contrast vector, it needs the names of the model matrix and everything should be 0 at default before adding the weights
  contrast_vec = rep(0, ncol(mod_mat))
  names(contrast_vec) = colnames(mod_mat)
  
  # combine the coefficient columns so we have the effects, then weight them because they're supposed to add to 1
  all_cols = c(first_col, other_cols)
  contrast_vec[all_cols] = 1 / length(all_cols)
  
  # extract results using numeric contrast vector, and get the shrunken L2FCs because they're more reliable/conservative
  res = results(INTWLD, contrast = contrast_vec, alpha = 0.05)
  res_shrk = lfcShrink(INTWLD, contrast = contrast_vec, type = "ashr", res = res)
  
  # put in lists because stageR needs stuff in lists
  p_list[[int]] = res$pvalue
  log2FC_list[[paste0(int, "_L2FC")]] = res$log2FoldChange
  log2FC_shrk_list[[paste0(int, "_shrk_L2FC")]] = res_shrk$log2FoldChange
}

# Restructuring to prepare for stage R
p_confirm = do.call(cbind, p_list)
rownames(p_confirm) = rownames(INTWLD)
l2fc_confirm = do.call(cbind, log2FC_list)
rownames(l2fc_confirm) = rownames(INTWLD)
l2fc_shrk_confirm = do.call(cbind, log2FC_shrk_list)
rownames(l2fc_shrk_confirm) = rownames(INTWLD)


# Now a stage R object will be created which will contain the information required for stageR analysis to be run in the format required (including the omnibus test results: pscreen, and the wald test results: pconfirm)
# I have also clarified to th stage R function that the p values for the LRT were NOT adjusted prior
stg_r_ob = stageR(
  pScreen = p_screen,
  pConfirmation = p_confirm,
  pScreenAdjusted = F
)

# Here we are performing the stagewise adjustment, this tests ACROSS all genes (FDR corrected with BH) for significance p.adj <0.05
# Genes which PASS this screening stage show significant synergy or antagonism associated with the interaction variables.
# Any genes that pass the screening stage are tested at the confirmation stage where FWER or FDR are controlled ACROSS hypotheses rather than genes this time
# This looks at each treatment comparison individually, and in this case we control the within-gene family wise error rate with Holm, alpha 0.05.
# This is a very stringent, conservative control of FWER, and controls the probability of making at least 1 false positive across all hypotheses.
# I'm going to do a version of this with BH later which will be used for pathways and exploratory analysis where we want to be less conservative, but the Holm method is great for when we want to look at specific genes/tests rather than pathways.
stg_r_ob = stageWiseAdjustment(stg_r_ob, method = "holm", alpha = 0.05)


# extract results
stage_r_res = getAdjustedPValues (stg_r_ob, order = F, onlySignificantGenes = F)

# Note for my output log
paste("Check the number of identified non-additive effects from any Interaction (Holm method):")
paste("WARNING: These are not necessarily statistically detectable in this analysis - that will depend on the LFC, which will be detected in the DEG_list_maker.R script")
sum(stage_r_res[ , -1] < 0.05, na.rm = TRUE)

# Note for output log
paste("Check number of genes where at least one Interaction displays non-additive effects (Holm method):")
paste("WARNING: These are not necessarily statistically detectable in this analysis - that will depend on the LFC, which will be detected in the DEG_list_maker.R script")
sum(apply(stage_r_res[ , -1], 1, function(row) any(row < 0.05, na.rm = TRUE)))

# Add the log 2 fold change values for the Wald tests to the results table
res4_tmp = merge(stage_r_res, l2fc_confirm, by = "row.names")
rownames(res4_tmp)=res4_tmp$Row.names
res4_tmp$Row.names = NULL
res4_table = merge(res4_tmp, l2fc_shrk_confirm, by = "row.names")
rownames(res4_table)=res4_table$Row.names
res4_table$Row.names = NULL

# Save results to a table
write.table(res4_table, "DEG_lists/res4_table.tsv", row.names= T, sep ="\t", quote = F)


## This is now going to be exactly the same thing but I will use BH - this isn't explicitly provided as a method WITHIN stageR, so I have to apply it to pconfirm outside of stage R and then select method="none".
# BH is less conservative it controls the PROPORTION of false positives and thus will be better for identifying DE pathways which could be more negatively affected by false negatives
BH_adjusted_matrix = t(apply(p_confirm, 1, function(p) p.adjust(p, method = "BH")))

stage_r_BH = stageR(pScreen = p_screen,
                    pConfirmation = BH_adjusted_matrix,
                    pScreenAdjusted = FALSE)

stage_r_BH = stageWiseAdjustment(object = stage_r_BH,
                                 method = "none", # Using none here because I have already BH adjusted in this case
                                 alpha = 0.05)

# Extract results
stage_r_BH_res = getAdjustedPValues (stage_r_BH, order = F, onlySignificantGenes = F)

# Note for output log
paste("Check the number of identified non-additive effects from any Interaction (BH method):")
paste("WARNING: These are not necessarily statistically detectable in this analysis - that will depend on the LFC, which will be detected in the DEG_list_maker.R script")
sum(stage_r_BH_res[ , -1] < 0.05, na.rm = TRUE)

# Note for output log
paste("Check number of genes where at least one Interaction displays non-additive effects (BH method):")
paste("WARNING: These are not necessarily statistically detectable in this analysis - that will depend on the LFC, which will be detected in the DEG_list_maker.R script")
sum(apply(stage_r_BH_res[ , -1], 1, function(row) any(row < 0.05, na.rm = TRUE)))

# Add the log 2 fold change values for the wald tests to the results table
res4_BH_tmp = merge(stage_r_BH_res, l2fc_confirm, by = "row.names")
rownames(res4_BH_tmp)=res4_BH_tmp$Row.names
res4_BH_tmp$Row.names = NULL
res4_BH_table = merge(res4_BH_tmp, l2fc_shrk_confirm, by = "row.names")
rownames(res4_BH_table)=res4_BH_table$Row.names
res4_BH_table$Row.names = NULL

# Save results
write.table(res4_BH_table, "DEG_lists/res4_BH_table.tsv", row.names= T, sep ="\t", quote = F)





#### 5. CHEMICAL INTERACTIONS PER TIME POINT ANALYSIS ####
# This version is for GSEA analysis over time, aka the middle ground: more detailed than the interaction overview but focuses on all interaction associated DE at each timepoint
# so uses LRT1 again to filter and not quite as zoomed in and sensitive to interaction jumps over time as the within interaction time comparisons are...

# Assign the UNadjusted p values from LRT1 to be the screening stage 1 for the stageR FDR correction, and provide the new object with the gene names
p_screen = resINTLRT1$pvalue
names(p_screen) = rownames(resINTLRT1)

# pull different timepoint types
timepoints = levels(INTWLD$Time)

# identify different interaction types
interactions = c("Azo_Pro", "Cyp_Pro", "Imid_Pro")

#lists for loop
p_list = list()
log2FC_list = list()
log2FC_shrk_list = list()

for (int in interactions) {
  
  # Extract the changing chem
  i = strsplit(int, "_Pro")[[1]]
  
  for (time in timepoints) {
    
    contrast_label = paste0(int, "_", time, "h")
    
    # build coefficient names exactly as in resultsNames(INTWLD)
    coef_1 = paste0(i, "pres.Propres.Time", time)
    coef_baseline = paste0(i, "pres.Propres")
    
    if (time == "1") {
      
      # When time is baseline, only use coef_baseline
      if (coef_baseline %in% resultsNames(INTWLD)) {
        res = results(INTWLD, name = coef_baseline, alpha = 0.05)
        res_shrk = lfcShrink(INTWLD, coef = coef_baseline, type="ashr", res=res)
        
        p_list[[contrast_label]] = res$pvalue
        log2FC_list[[paste0(contrast_label, "_L2FC")]] = res$log2FoldChange
        log2FC_shrk_list[[paste0(contrast_label, "_shrk_L2FC")]] = res_shrk$log2FoldChange
        
      } else {
        message(paste("Skipping", contrast_label, "- coefficient not found"))
      }
    } else {
      # time point are not baseline: contrast of two coefficients
      
      if (all(c(coef_1, coef_baseline) %in% resultsNames(INTWLD))) {
        res = results(INTWLD, contrast = list(c(coef_baseline, coef_1)), alpha = 0.05)
        res_shrk = lfcShrink(INTWLD, contrast = list(c(coef_baseline, coef_1)), type="ashr", res=res)
        
        p_list[[contrast_label]] = res$pvalue
        log2FC_list[[paste0(contrast_label, "_L2FC")]] = res$log2FoldChange
        log2FC_shrk_list[[paste0(contrast_label, "_shrk_L2FC")]] = res_shrk$log2FoldChange
      } else {
        message(paste("Skipping", contrast_label, "- coefficients not found"))
      }
    }
  }
}


p_confirm = do.call(cbind, p_list)
rownames(p_confirm) = rownames(INTWLD)
l2fc_confirm = do.call(cbind, log2FC_list)
rownames(l2fc_confirm) = rownames(INTWLD)
l2fc_shrk_confirm = do.call(cbind, log2FC_shrk_list)
rownames(l2fc_shrk_confirm) = rownames(INTWLD)

stg_r_ob = stageR(
  pScreen = p_screen,
  pConfirmation = p_confirm,
  pScreenAdjusted = F
)

stg_r_ob = stageWiseAdjustment(stg_r_ob, method = "holm", alpha = 0.05)


# Extract results
stage_r_res = getAdjustedPValues (stg_r_ob, order = F, onlySignificantGenes = F)

# Note for output log
paste("Check the number of identified non-additive effects from any Interaction at any timepoint (Holm method):")
paste("WARNING: These are not necessarily statistically detectable in this analysis - that will depend on the LFC, which will be detected in the DEG_list_maker.R script")
sum(stage_r_res[ , -1] < 0.05, na.rm = TRUE)

# Note for output log
paste("Check number of genes where at least one Interaction at any timepoint displays non-additive effects (BH method):")
paste("WARNING: These are not necessarily statistically detectable in this analysis - that will depend on the LFC, which will be detected in the DEG_list_maker.R script")
sum(apply(stage_r_res[ , -1], 1, function(row) any(row < 0.05, na.rm = TRUE)))

# Add the log 2 fold change values for the wald tests to the results table
res5_tmp = merge(stage_r_res, l2fc_confirm, by = "row.names")
rownames(res5_tmp)=res5_tmp$Row.names
res5_tmp$Row.names = NULL
res5_table = merge(res5_tmp, l2fc_shrk_confirm, by = "row.names")
rownames(res5_table)=res5_table$Row.names
res5_table$Row.names = NULL

# Save results
write.table(res5_table, "DEG_lists/res5_table.tsv", row.names= T, sep ="\t", quote = F)


## BH edit
BH_adjusted_matrix = t(apply(p_confirm, 1, function(p) p.adjust(p, method = "BH")))

stage_r_BH = stageR(pScreen = p_screen,
                    pConfirmation = BH_adjusted_matrix,
                    pScreenAdjusted = FALSE)

stage_r_BH = stageWiseAdjustment(object = stage_r_BH,
                                 method = "none", # Using none here because I have already BH adjusted in this case
                                 alpha = 0.05)

# Extract results
stage_r_BH_res = getAdjustedPValues (stage_r_BH, order = F, onlySignificantGenes = F)

# Note for output log
paste("Check the number of identified non-additive effects from any Interaction at any timepoint (BH method):")
paste("WARNING: These are not necessarily statistically detectable in this analysis - that will depend on the LFC, which will be detected in the DEG_list_maker.R script")
sum(stage_r_BH_res[ , -1] < 0.05, na.rm = TRUE)

# Note for output log
paste("Check number of genes where at least one Interaction at any timepoint displays non-additive effects (BH method):")
paste("WARNING: These are not necessarily statistically detectable in this analysis - that will depend on the LFC, which will be detected in the DEG_list_maker.R script")
sum(apply(stage_r_BH_res[ , -1], 1, function(row) any(row < 0.05, na.rm = TRUE)))

# Add the log 2 fold change values for the wald tests to the results table
res5_BH_tmp = merge(stage_r_BH_res, l2fc_confirm, by = "row.names")
rownames(res5_BH_tmp)=res5_BH_tmp$Row.names
res5_BH_tmp$Row.names = NULL
res5_BH_table = merge(res5_BH_tmp, l2fc_shrk_confirm, by = "row.names")
rownames(res5_BH_table)=res5_BH_table$Row.names
res5_BH_table$Row.names = NULL

# Save results
write.table(res5_BH_table, "DEG_lists/res5_BH_table.tsv", row.names= T, sep ="\t", quote = F)





#### 6. TIME DYNAMICS WITHIN CHEMICAL INTERACTIONS ANALYSIS ####

# Here what I'm doing is very similar to the previous method but here genes which pass this screening stage show significant interaction effects associated with the interaction:time interaction term.
# This means that they show both synergy and antagonism CHANGES due to time-dependent interaction affects, which will identify interesting interaction based behaviour of gene expression over time.
# We will then be able to compare individual time points WITHIN each interaction.
p_screen = resINTLRT2$pvalue
names(p_screen) = rownames(resINTLRT2)

# Loop through all of the INTERACTION Wald tests to be performed and extract the p values and log fold changes for stage 2 FWER correction.

# identify different interaction types
interactions = c("Azo_Pro", "Cyp_Pro", "Imid_Pro")

timepoints = levels(INTWLD$Time)
time_pairs = combn(timepoints, 2, simplify = FALSE)

p_list = list()
log2FC_list = list()
log2FC_shrk_list = list()

for (int in interactions) {
  
  # Extract the changing chem
  i = strsplit(int, "_Pro")[[1]]

  for (pair in time_pairs) {
    time1 = pair[1]
    time2 = pair[2]
    contrast_label = paste0(int, "_at_", time2, "vs", time1)
    
    # interaction coefficients
    coef1 = paste0(i, "pres.Propres.Time", time2)  # interaction at time2
    coef2 = paste0(i, "pres.Propres.Time", time1)  # interaction at time1
    
    if (time1 == timepoints[1]) {   # baseline = first level, e.g. "1"
      # only need to do coef 1
      if (all(c(coef1) %in% resultsNames(INTWLD))) {
        res = results(INTWLD, contrast = list(c(coef1)), alpha = 0.05)
        res_shrk = lfcShrink(INTWLD, contrast = list(c(coef1)), type="ashr", res=res)
        p_list[[contrast_label]] = res$pvalue
        log2FC_list[[paste0(contrast_label, "_L2FC")]] = res$log2FoldChange
        log2FC_shrk_list[[paste0(contrast_label, "_shrk_L2FC")]] = res_shrk$log2FoldChange
      } else {
        message("Skipping ", contrast_label, " - needed coefficients missing")
      }
    } else {
      # General case: (Time_t2 + T:Time_t2) - (Time_t1 + T:Time_t1)
      if (all(c(coef1, coef2) %in% resultsNames(INTWLD))) {
        res = results(INTWLD, contrast = list(c(coef1), c(coef2)), alpha = 0.05)
        res_shrk = lfcShrink(INTWLD, contrast = list(c(coef1), c(coef2)), type="ashr", res=res)
        p_list[[contrast_label]] = res$pvalue
        log2FC_list[[paste0(contrast_label, "_L2FC")]] = res$log2FoldChange
        log2FC_shrk_list[[paste0(contrast_label, "_shrk_L2FC")]] = res_shrk$log2FoldChange
      } else {
        message("Skipping ", contrast_label, " - needed coefficients missing")
      }
    }
  }
}


p_confirm = do.call(cbind, p_list)
rownames(p_confirm) = rownames(INTWLD)
l2fc_confirm = do.call(cbind, log2FC_list)
rownames(l2fc_confirm) = rownames(INTWLD)
l2fc_shrk_confirm = do.call(cbind, log2FC_shrk_list)
rownames(l2fc_shrk_confirm) = rownames(INTWLD)

stg_r_ob = stageR(
  pScreen = p_screen,
  pConfirmation = p_confirm,
  pScreenAdjusted = F
)

stg_r_ob = stageWiseAdjustment(stg_r_ob, method = "holm", alpha = 0.05)


# Extract results
stage_r_res = getAdjustedPValues (stg_r_ob, order = F, onlySignificantGenes = F)

# Note for output log
paste("Check the number of instances where the interaction effect CHANGES between timepoints (aka becomes more or less synergytic/antagonistic) (Holm method):")
paste("WARNING: These are not necessarily statistically detectable in this analysis - that will depend on the LFC, which will be detected in the DEG_list_maker.R script")
sum(stage_r_res[ , -1] < 0.05, na.rm = TRUE)

# Note for output log
paste("Check the number of genes where at least one interaction effect CHANGES between timepoints (aka becomes more or less synergytic/antagonistic) (Holm method):")
paste("WARNING: These are not necessarily statistically detectable in this analysis - that will depend on the LFC, which will be detected in the DEG_list_maker.R script")
sum(apply(stage_r_res[ , -1], 1, function(row) any(row < 0.05, na.rm = TRUE)))

# Add the log 2 fold change values for the wald tests to the results table
res6_tmp = merge(stage_r_res, l2fc_confirm, by = "row.names")
rownames(res6_tmp)=res6_tmp$Row.names
res6_tmp$Row.names = NULL
res6_table = merge(res6_tmp, l2fc_shrk_confirm, by = "row.names")
rownames(res6_table)=res6_table$Row.names
res6_table$Row.names = NULL

# Save results
write.table(res6_table, "DEG_lists/res6_table.tsv", row.names= T, sep ="\t", quote = F)


## BH edit
BH_adjusted_matrix = t(apply(p_confirm, 1, function(p) p.adjust(p, method = "BH")))

stage_r_BH = stageR(pScreen = p_screen,
                    pConfirmation = BH_adjusted_matrix,
                    pScreenAdjusted = FALSE)

stage_r_BH = stageWiseAdjustment(object = stage_r_BH,
                                 method = "none", # Using none here because I have already BH adjusted in this case
                                 alpha = 0.05)

# Extract results
stage_r_BH_res = getAdjustedPValues (stage_r_BH, order = F, onlySignificantGenes = F)

# Note for output log
paste("Check the number of instances where the interaction effect CHANGES between timepoints (aka becomes more or less synergytic/antagonistic) (BH method):")
paste("WARNING: These are not necessarily statistically detectable in this analysis - that will depend on the LFC, which will be detected in the DEG_list_maker.R script")
sum(stage_r_BH_res[ , -1] < 0.05, na.rm = TRUE)

# Note for output log
paste("Check the number of genes where at least one interaction effect CHANGES between timepoints (aka becomes more or less synergytic/antagonistic)  (BH method):")
paste("WARNING: These are not necessarily statistically detectable in this analysis - that will depend on the LFC, which will be detected in the DEG_list_maker.R script")
sum(apply(stage_r_BH_res[ , -1], 1, function(row) any(row < 0.05, na.rm = TRUE)))

# Add the log 2 fold change values for the wald tests to the results table
res6_BH_tmp = merge(stage_r_BH_res, l2fc_confirm, by = "row.names")
rownames(res6_BH_tmp)=res6_BH_tmp$Row.names
res6_BH_tmp$Row.names = NULL
res6_BH_table = merge(res6_BH_tmp, l2fc_shrk_confirm, by = "row.names")
rownames(res6_BH_table)=res6_BH_table$Row.names
res6_BH_table$Row.names = NULL

# Save results
write.table(res6_BH_table, "DEG_lists/res6_BH_table.tsv", row.names= T, sep ="\t", quote = F)





#### Save ####

# Save .RData file
save.image(file=paste0("Filtered_DESeq2_Results", ".RData"))