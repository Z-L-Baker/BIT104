#### R SCRIPT FOR MAKING VOLCANO PLOTS FROM DESEQ2 OUTPUT DATA ####

#### Set-up ####
rm(list=ls())
#Set working directory
workDir = "/mnt/clusters/admiral/data/c21082179/BIT104/New_Analysis"  # working directory for the R session
setwd(workDir)

## Read in the .RData file from the omnibus_method_results.R script :>
load("Filtered_DESeq2_Results.RData")

# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# 
# BiocManager::install('EnhancedVolcano')
# devtools::install_github('kevinblighe/EnhancedVolcano')

prot_names = read.table("metadata/gene2protID.tsv", header=F, sep="\t", na.strings="")
rownames(prot_names) = prot_names[,1]
prots = prot_names[ -c(1,3:4) ]
colnames(prots) = "protID"
                    
library(EnhancedVolcano)

dir.create(paste0("figures/volcano_plots"))

VULCAN = function(df, l2f_min, subfolder){
  
  # Make the directories needed
  #dir.create(paste0("figures/volcano_plots/", subfolder))
  dir.create(paste0("figures/volcano_plots/", subfolder, "_shrknLFC"))
  
  # take columns not ending in L2FC (aka take all the pval columns and none of the log 2 fold change columns shrunken or otherwise)
  pval_cols = grep("_.*$", names(df), value = TRUE)
  pval_cols = pval_cols[!grepl("_L2FC$", pval_cols)]  # exclude L2FC columns
  
  for (i in seq_along(pval_cols)) {
    pval_col = pval_cols[i]
    l2fc_col = paste0(pval_col, "_shrk", "_L2FC")
    
    if (l2fc_col %in% names(df)) {
      pvals = df[[pval_col]]
      l2fc  = df[[l2fc_col]]
      
      tmp_df = data.frame(
        log2FC = l2fc,
        pval = pvals,
        row.names = rownames(df)
      )
      
      tmp_df = merge(tmp_df, prots, by = "row.names")
      rownames(tmp_df) = tmp_df$Row.names 
      tmp_df = tmp_df[ -c(1) ]
      
      tmp_df = na.omit(tmp_df)
      
      # make a ranking score because too many genes are getting labelled and it's messy to say the least
      pval_cutoff = 0.05
      fc_cutoff   = l2f_min
      
      tmp_df$rankScore = ifelse(
        !is.na(tmp_df$pval) & !is.na(tmp_df$log2FC) &
          abs(tmp_df$log2FC) > fc_cutoff &
          tmp_df$pval < pval_cutoff,
        sqrt( (abs(tmp_df$log2FC))^2 + (-log10(tmp_df$pval))^2 ),
        0
      )
      
      # only positive scores (for cases where nothing is actually passing the thresholds)
      tmp_df = tmp_df[order(tmp_df$rankScore, decreasing = TRUE), ]
      topGenes = head(tmp_df$protID, 5)
      
      p = EnhancedVolcano(tmp_df,
                          lab = ifelse(tmp_df$protID %in% topGenes, tmp_df$protID, ""),  # only label topGenes
                          x = 'log2FC',
                          y = 'pval',
                          #xlim = c(-7.5, 7.5),
                          #ylim = c(0, 70),
                          axisLabSize = 20,
                          subtitle = pval_col,
                          pCutoff = 0.05,
                          FCcutoff = l2f_min,
                          cutoffLineCol = 'grey30',
                          pointSize = 3.5,
                          labSize = 5.0,
                          col=c('grey50', 'grey50', 'grey50', 'red3'),
                          drawConnectors = TRUE,
                          widthConnectors = 0.5,
                          colConnectors = 'black',
                          arrowheads = FALSE,
                          colAlpha = 0.5,
                          legendPosition = 'none',
                          selectLab = NULL
      )
      
      ggsave(paste0("figures/volcano_plots/", subfolder, "_shrknLFC/", pval_col, ".png"), plot = p, width = 12, height = 10)
    }
    
    #l2fc_col = paste0(pval_col, "_L2FC")
    
    #if (l2fc_col %in% names(df)) {
    #  pvals = df[[pval_col]]
    #  l2fc  = df[[l2fc_col]]
    #  
    #  tmp_df = data.frame(
    #    log2FC = l2fc,
    #    pval = pvals,
    #    row.names = rownames(df)
    #  )
    #  
    #  tmp_df = merge(tmp_df, prots, by = "row.names")
    #  rownames(tmp_df) = tmp_df$Row.names 
    #  tmp_df = tmp_df[ -c(1) ]
    #  
    #  tmp_df = na.omit(tmp_df)
    #  
    #  # make a ranking score because too many genes are getting labelled and it's messy to say the least
    #  pval_cutoff = 0.05
    #  fc_cutoff   = l2f_min
    #  
    #  tmp_df$rankScore = ifelse(
    #    !is.na(tmp_df$pval) & !is.na(tmp_df$log2FC) &
    #      abs(tmp_df$log2FC) > fc_cutoff &
    #      tmp_df$pval < pval_cutoff,
    #    sqrt( (abs(tmp_df$log2FC))^2 + (-log10(tmp_df$pval))^2 ),
    #    0
    #  )
    #  
    #  # only positive scores (for cases where nothing is actually passing the thresholds)
    #  tmp_df = tmp_df[order(tmp_df$rankScore, decreasing = TRUE), ]
    #  topGenes = head(tmp_df$protID, 5)
    #  
    #  p = EnhancedVolcano(tmp_df,
    #                      lab = ifelse(tmp_df$protID %in% topGenes, tmp_df$protID, ""),  # only label topGenes
    #                      x = 'log2FC',
    #                      y = 'pval',
    #                      #xlim = c(-7.5, 7.5),
    #                      #ylim = c(0, 70),
    #                      subtitle = pval_col,
    #                      pCutoff = 0.05,
    #                      FCcutoff = l2f_min,
    #                      cutoffLineCol = 'grey30',
    #                      pointSize = 3.0,
    #                      labSize = 4.0,
    #                      col=c('grey50', 'grey50', 'grey50', 'red3'),
    #                      drawConnectors = TRUE,
    #                      widthConnectors = 0.5,
    #                      colConnectors = 'black',
    #                      arrowheads = FALSE,
    #                      colAlpha = 0.5,
    #                      legendPosition = 'none',
    #                      selectLab = NULL
    #  )
    #  
    #  ggsave(paste0("figures/volcano_plots/", subfolder, "/", pval_col, ".png"), plot = p, width = 12, height = 10)
    #}
  }
}

VULCAN(df = res1_table, l2f_min = 1, subfolder = "TreatmentAverage_Holm_2fold")
VULCAN(res1_BH_table, 1, "TreatmentAverage_BH_2fold")
VULCAN(res2_table, 2, "TreatmentPerTime_Holm_4fold")
VULCAN(res2_BH_table, 2, "TreatmentPerTime_BH_4fold")
VULCAN(res3_table, 2, "TimeWithinTreatment_Holm_4fold")
VULCAN(res3_BH_table, 2, "TimeWithinTreatment_BH_4fold")

VULCAN(df = res4_table, l2f_min = 1, subfolder = "InteractionAverage_Holm_2fold")
VULCAN(res4_BH_table, 1, "InteractionAverage_BH_2fold")
VULCAN(res5_table, 2, "InteractionPerTime_Holm_4fold")
VULCAN(res5_BH_table, 2, "InteractionPerTime_BH_4fold")
VULCAN(res6_table, 2, "TimeWithinInteraction_Holm_4fold")
VULCAN(res6_BH_table, 2, "TimeWithinInteraction_BH_4fold")

##### GIFS ####
#
#dir.create("figures/volcano_plots/gifs")
#
#
### Azo-Pro###
#png_files = list.files(
#  "figures/volcano_plots/TreatmentPerTime_Holm_4fold_shrknLFC",
#  pattern = "^Azo-Pro.*\\.png$",
#  full.names = TRUE
#)
#nums = as.numeric(gsub("^Azo-Pro_([0-9]+).*\\.png$", "\\1", basename(png_files)))
#png_files = png_files[order(nums)]
#
#gifski(png_files, gif_file = "figures/volcano_plots/gifs/Azo-Pro.gif", delay = 1.5)  # delay in seconds per frame
#
#
### Azoxy ###
#png_files = list.files(
#  "figures/volcano_plots/TreatmentPerTime_Holm_4fold_shrknLFC",
#  pattern = "^Azoxy.*\\.png$",
#  full.names = TRUE
#)
#nums = as.numeric(gsub("^Azoxy_([0-9]+).*\\.png$", "\\1", basename(png_files)))
#png_files = png_files[order(nums)]
#
#gifski(png_files, gif_file = "figures/volcano_plots/gifs/Azoxy.gif", delay = 1.5)  # delay in seconds per frame
#
### Cyp ###
#png_files = list.files(
#  "figures/volcano_plots/TreatmentPerTime_Holm_4fold_shrknLFC",
#  pattern = "^Cyp_.*\\.png$",
#  full.names = TRUE
#)
#nums = as.numeric(gsub("^Cyp_([0-9]+).*\\.png$", "\\1", basename(png_files)))
#png_files = png_files[order(nums)]
#
#gifski(png_files, gif_file = "figures/volcano_plots/gifs/Cyp.gif", delay = 1.5)  # delay in seconds per frame
#
### CypPro ###
#png_files = list.files(
#  "figures/volcano_plots/TreatmentPerTime_Holm_4fold_shrknLFC",
#  pattern = "^Cyp-Pro.*\\.png$",
#  full.names = TRUE
#)
#nums = as.numeric(gsub("^Cyp-Pro_([0-9]+).*\\.png$", "\\1", basename(png_files)))
#png_files = png_files[order(nums)]
#
#gifski(png_files, gif_file = "figures/volcano_plots/gifs/Cyp-Pro.gif", delay = 1.5)  # delay in seconds per frame
#
### ImdPro ###
#png_files = list.files(
#  "figures/volcano_plots/TreatmentPerTime_Holm_4fold_shrknLFC",
#  pattern = "^Imd-Pro.*\\.png$",
#  full.names = TRUE
#)
#nums = as.numeric(gsub("^Imd-Pro_([0-9]+).*\\.png$", "\\1", basename(png_files)))
#png_files = png_files[order(nums)]
#
#gifski(png_files, gif_file = "figures/volcano_plots/gifs/Imd-Pro.gif", delay = 1.5)  # delay in seconds per frame
#
### Imid ###
#png_files = list.files(
#  "figures/volcano_plots/TreatmentPerTime_Holm_4fold_shrknLFC",
#  pattern = "^Imid_.*\\.png$",
#  full.names = TRUE
#)
#nums = as.numeric(gsub("^Imid_([0-9]+).*\\.png$", "\\1", basename(png_files)))
#png_files = png_files[order(nums)]
#
#gifski(png_files, gif_file = "figures/volcano_plots/gifs/Imid.gif", delay = 1.5)  # delay in seconds per frame
#
### Pro ###
#png_files = list.files(
#  "figures/volcano_plots/TreatmentPerTime_Holm_4fold_shrknLFC",
#  pattern = "^Pro_.*\\.png$",
#  full.names = TRUE
#)
#nums = as.numeric(gsub("^Pro_([0-9]+).*\\.png$", "\\1", basename(png_files)))
#png_files = png_files[order(nums)]
#
#gifski(png_files, gif_file = "figures/volcano_plots/gifs/Pro.gif", delay = 1.5)  # delay in seconds per frame
#
