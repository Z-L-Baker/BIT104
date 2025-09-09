#### Set-up ####
rm(list=ls())
#Set working directory
workDir = "/mnt/clusters/admiral/data/c21082179/BIT104/New_Analysis"  # working directory for the R session
setwd(workDir)



## Read in the .RData file from the omnibus_method_results.R script :>
load("Filtered_DESeq2_Results.RData")

dir.create("DEG_lists/")

#### I'm going toput all of this in a function because I don't want to nest too many loops

DEG_LIST_MAKER = function(df, l2f_min, subfolder, shrk = TRUE) {
  
  if (shrk == TRUE){
    print("Function will run with shrunken L2FC values")
    shrknLFC = "_shrk"
  } else {
    print("Function will run with non-shrunken L2FC values")
    shrknLFC = ""
  }
  
  # Make the directories needed
  dir.create(paste0("DEG_lists/", subfolder, shrknLFC))
  dir.create(paste0("DEG_lists/", subfolder, shrknLFC, "/all_DEGs/"))
  dir.create(paste0("DEG_lists/", subfolder, shrknLFC, "/upreg_only/"))
  dir.create(paste0("DEG_lists/", subfolder, shrknLFC, "/downreg_only/"))
  
  # take columns not ending in L2FC (aka take all the pval columns and none of the log 2 fold change columns shrunken or otherwise)
  pval_cols = grep("_.*$", names(df), value = TRUE)
  pval_cols = pval_cols[!grepl("_L2FC$", pval_cols)]  # exclude L2FC columns
  
  # make empty matrix which will later contain true or false condition based on if the L2FC is wihtin the set threshold (which as prior determined will be 2 for time treatment and 1 for treatment)
  # it'll be the same size as the results dataframe and set to false by default, only resulting in a true if the L2FC and p value for that hypothesis and gene meet the thresholds :>
  condition_matrix = matrix(FALSE, nrow = nrow(df), ncol = length(pval_cols))
  colnames(condition_matrix) = pval_cols
  rownames(condition_matrix) = rownames(df)
  
  ## also make matrices for SPECIFIC upreg and downreg data
  condition_matrix_up = matrix(FALSE, nrow = nrow(df), ncol = length(pval_cols))
  colnames(condition_matrix_up) = pval_cols
  rownames(condition_matrix_up) = rownames(df)
  
  condition_matrix_down = matrix(FALSE, nrow = nrow(df), ncol = length(pval_cols))
  colnames(condition_matrix_down) = pval_cols
  rownames(condition_matrix_down) = rownames(df)
  
  
  # iterate through genes and fill TRUE where the threshold condition is met
  for (i in seq_along(pval_cols)) {
    pval_col = pval_cols[i]
    l2fc_col = paste0(pval_col, shrknLFC, "_L2FC")
    
    if (l2fc_col %in% names(df)) {
      pvals = df[[pval_col]]
      l2fc  = df[[l2fc_col]]
      
      condition_matrix[, i] = pvals < 0.05 & abs(l2fc) > l2f_min
    }
  }
  
  # do the same iteration technique but for upregulated genes specifically so we can make lists from them
  for (i in seq_along(pval_cols)) {
    pval_col = pval_cols[i]
    l2fc_col = paste0(pval_col, shrknLFC, "_L2FC")
    
    if (l2fc_col %in% names(df)) {
      pvals = df[[pval_col]]
      l2fc  = df[[l2fc_col]]
      
      condition_matrix_up[, i] = pvals < 0.05 & l2fc > l2f_min
    }
  }
  
  # an for downregulated genes specifically 
  for (i in seq_along(pval_cols)) {
    pval_col = pval_cols[i]
    l2fc_col = paste0(pval_col, shrknLFC, "_L2FC")
    
    if (l2fc_col %in% names(df)) {
      pvals = df[[pval_col]]
      l2fc  = df[[l2fc_col]]
      
      condition_matrix_down[, i] = pvals < 0.05 & l2fc < -(l2f_min)
    }
  }
  
  
  # extract total number of true instances (number of DE instances above pval and L2FC threshold)
  total_instances = sum(condition_matrix, na.rm = TRUE)
  
  # extract number of rows with at least one true instance (number of GENES which show DE in at least one instance between treatment conditions)
  rows_with_instance = sum(apply(condition_matrix, 1, any, na.rm = TRUE))
  
  # extract number of genes with at least one true instance (instance of upreg DE or downreg DE) between any hypothesis treatment comparison
  overall_significant_genes = rownames(condition_matrix)[apply(condition_matrix, 1, any, na.rm = TRUE)]
  overall_significant_up_genes = rownames(condition_matrix_up)[apply(condition_matrix_up, 1, any, na.rm = TRUE)]
  overall_significant_down_genes = rownames(condition_matrix_down)[apply(condition_matrix_down, 1, any, na.rm = TRUE)]
  
  # most importantly: extract gene names per treatment comparison and save these toa list (will later save as DEG list text files)
  comparison_gene_lists = lapply(pval_cols, function(col) {
    rownames(condition_matrix)[which(condition_matrix[, col])]
  })
  names(comparison_gene_lists) = pval_cols
  
  comparison_up_gene_lists = lapply(pval_cols, function(col) {
    rownames(condition_matrix_up)[which(condition_matrix_up[, col])]
  })
  names(comparison_up_gene_lists) = pval_cols
  
  comparison_down_gene_lists = lapply(pval_cols, function(col) {
    rownames(condition_matrix_down)[which(condition_matrix_down[, col])]
  })
  names(comparison_down_gene_lists) = pval_cols
  
  
  # REPORT STATS FOR LOG FILE
  cat("Number of DE instances:", total_instances, "\n")
  cat("Number of Genes which show DE in at least one instance:", rows_with_instance, "\n")
  
  # Not needed now I've merged them
  #pvals = df[["Carrier_vs_Water"]]
  #l2fc  = df[["Carrier_vs_Water_L2FC"]]
  #valid_rows = !is.na(pvals) & pvals < 0.05 & abs(l2fc) > 1
  #count = sum(valid_rows, na.rm = TRUE)
  #cat("Number of Genes which show DE between Carrier and Water treatment groups:", count, "\n")
  
  
  # Optional: view counts
  cat("Number of upregulated genes with at least one significant comparison:", length(overall_significant_up_genes), "\n")
  cat("Genes per comparison:\n")
  print(sapply(comparison_up_gene_lists, length))
  
  cat("Number of downregulated genes with at least one significant comparison:", length(overall_significant_down_genes), "\n")
  cat("Genes per comparison:\n")
  print(sapply(comparison_down_gene_lists, length))
  
  
  # FINALLY SAVE THE DEG LISTS TO FILES FOR DOWNSTREAM ANALSYSES :>
  for (i in 1:length(comparison_gene_lists)){
    genes = comparison_gene_lists[[i]]
    name = names(comparison_gene_lists)[i]
    write(genes, file = paste0("DEG_lists/", subfolder, shrknLFC, "/all_DEGs/", name, ".txt"), sep =" ")
  }
  
  for (i in 1:length(comparison_up_gene_lists)){
    genes = comparison_up_gene_lists[[i]]
    name = names(comparison_up_gene_lists)[i]
    write(genes, file = paste0("DEG_lists/", subfolder, shrknLFC, "/upreg_only/", name, ".txt"), sep =" ")
  }
  
  for (i in 1:length(comparison_down_gene_lists)){
    genes = comparison_down_gene_lists[[i]]
    name = names(comparison_down_gene_lists)[i]
    write(genes, file = paste0("DEG_lists/", subfolder, shrknLFC, "/downreg_only/", name, ".txt"), sep =" ")
  }
}


#### run the function on each dataframe, setting the l2fc as appropriate for that table and being transparent on FDR/FWER method (and l2fc for time treatment)

# I also want results for the 5-fold effect for time treatment just in case it actually does seem to give better results than 4-fold
# Although I still think 4-fold is justifiable and will allow for many more DE cases detected, especially considering even with 5fold we still have a couple of underpowered tests :<
# regardless, reporting both would be better form

# also if not wanting shrunken LFCs, please adapt the below to include shrk=F

cat("\nWorking on Treatment Time-Average Data (Shrunken LFCs):\n")
DEG_LIST_MAKER(df = res1_table, l2f_min = 1, subfolder = "TreatmentAverage_Holm_2fold")
DEG_LIST_MAKER(res1_BH_table, 1, "TreatmentAverage_BH_2fold")

cat("\nWorking on Treatment per Time Data (Shrunken LFCs):\n")
DEG_LIST_MAKER(res2_table, 2, "TreatmentPerTime_Holm_4fold")
DEG_LIST_MAKER(res2_BH_table, 2, "TreatmentPerTime_BH_4fold")

cat("\nWorking on Within-Treatment Time Data (Shrunken LFCs):\n")
DEG_LIST_MAKER(res3_table, 2, "TimeWithinTreatment_Holm_4fold")
DEG_LIST_MAKER(res3_BH_table, 2, "TimeWithinTreatment_BH_4fold")


cat("\nWorking on Interaction Time-Average Data (Shrunken LFCs):\n")
DEG_LIST_MAKER(df = res4_table, l2f_min = 1, subfolder = "InteractionAverage_Holm_2fold")
DEG_LIST_MAKER(res4_BH_table, 1, "InteractionAverage_BH_2fold")

cat("\nWorking on Interaction per Time Data (Shrunken LFCs):\n")
DEG_LIST_MAKER(res5_table, 2, "InteractionPerTime_Holm_4fold")
DEG_LIST_MAKER(res5_BH_table, 2, "InteractionPerTime_BH_4fold")

cat("\nWorking on Within-Interaction Time Data (Shrunken LFCs):\n")
DEG_LIST_MAKER(res6_table, 2, "TimeWithinInteraction_Holm_4fold")
DEG_LIST_MAKER(res6_BH_table, 2, "TimeWithinInteraction_BH_4fold")


#cat("\nWorking on Treatment Time-Average Data (Non-Shrunken LFCs):\n")
#DEG_LIST_MAKER(df = res1_table, l2f_min = 1, subfolder = "TreatmentAverage_Holm_2fold", shrk = F)
#DEG_LIST_MAKER(res1_BH_table, 1, "TreatmentAverage_BH_2fold", shrk = F)

#cat("\nWorking on Within-Treatment Time Data (Non-Shrunken LFCs):\n")
#DEG_LIST_MAKER(res2_table, 2, "TimeWithinTreatment_Holm_4fold", shrk = F)
#DEG_LIST_MAKER(res2_BH_table, 2, "TimeWithinTreatment_BH_4fold", shrk = F)

#cat("\nWorking on Treatment per Time Data (Non-Shrunken LFCs):\n")
#DEG_LIST_MAKER(res3_table, 2, "TreatmentPerTime_Holm_4fold", shrk = F)
#DEG_LIST_MAKER(res3_BH_table, 2, "TreatmentPerTime_BH_4fold", shrk = F)