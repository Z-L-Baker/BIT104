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
## BiocManager::install('ggridges')


library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr)
library(pheatmap)
library(tibble)

prot_names = read.table("metadata/gene2protID.tsv", header=F, sep="\t", na.strings="")
rownames(prot_names) = prot_names[,1]
prots = prot_names[ -c(1,3:4) ]
colnames(prots) = "protID"

dir.create(paste0("figures/volcano_ridge_plots"))
dir.create(paste0("figures/volcano_heatmap_plots"))

df = res2_table
for(trt in treatment_coefs){
  df_filtered = df %>%
    select(matches(paste0("^",trt)))

  df_final = df_filtered %>% select(c(1:6,13:18))
  df_final = merge(df_final, prots, by = "row.names")
  df_final = df_final %>%
    select(!c(Row.names))

  colnames(df_final) = c("T1_pval", "T4_pval", "T8_pval", "T12_pval", "T16_pval", "T24_pval", "T1_log2FC" , "T4_log2FC" , "T8_log2FC", "T12_log2FC", "T16_log2FC", "T24_log2FC", "gene")

  df_long = df_final %>%
    pivot_longer(
      cols = -gene,
      names_to = c("timepoint", ".value"),
      names_pattern = "(T\\d+)_(.*)"        # matches T1_log2FC → timepoint = T1, value = log2FC
    ) %>%
    mutate(negLogP = -log10(pval))

  df_long = na.omit(df_long)
  df_long$timepoint = factor(df_long$timepoint,
                            levels = c("T1", "T2", "T4", "T8", "T12", "T16", "T24"))

  p1 = ggplot(df_long, aes(x = log2FC, y = timepoint, color = negLogP, size = negLogP)) +
    ggtitle(paste0(trt)) +
    geom_point(alpha = 0.6) +
    geom_vline(xintercept = c(-2, 2), color = "black", linetype = "dashed") +
    scale_color_viridis_c(option = "C") +
    theme_minimal() +
    labs(x = "log2 Fold Change", y = "Timepoint (hours)",
        color = "-log10 p", size = "-log10 p")

    ggsave(paste0("figures/volcano_ridge_plots/", trt, "_Holm_4fold.png"), p1, width = 12, height = 5)

  
  sig_thresh = 0.05
  lfc_thresh = 2
  top_n = 20
  num_timepoints = length(unique(df_long$timepoint))
  
  # Handle duplicates
  df_long_unique = df_long %>%
    group_by(gene) %>%
    mutate(gene_count = n() / num_timepoints,
           like_index = ifelse(gene_count > 1,
                               ceiling(row_number() / num_timepoints),
                               NA),
           gene = ifelse(gene_count > 1,
                         paste0(gene, "-like", like_index),
                         gene)) %>%
    ungroup() %>%
    select(-gene_count, -like_index) %>%
    mutate(timepoint = factor(timepoint,
                              levels = c("T1", "T2", "T4", "T8", "T12", "T16", "T24")))
  
  # Identify top 10 per timepoint by pval < 0.05 and LFC > threshold
  top_genes_df = df_long_unique %>%
    filter(pval < sig_thresh & abs(log2FC) > lfc_thresh) %>%
    group_by(timepoint) %>%
    arrange(pval, desc(abs(log2FC))) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    mutate(top10 = TRUE)
  
  top_genes = unique(top_genes_df$gene)
  
  # Prepare heatmap matrix
  mat_heat = df_long_unique %>%
    filter(gene %in% top_genes) %>%
    select(gene, timepoint, log2FC) %>%
    pivot_wider(names_from = timepoint, values_from = log2FC) %>%
    column_to_rownames("gene") %>%
    as.matrix()
  
  # Prepare overlay (X for non-significant, * for significant LFC>2)
  overlay_mat = df_long_unique %>%
    filter(gene %in% top_genes) %>%
    mutate(mark = case_when(
      pval > sig_thresh ~ "X",
      abs(log2FC) > lfc_thresh & pval <= sig_thresh ~ "*",
      TRUE ~ ""
    )) %>%
    select(gene, timepoint, mark) %>%
    pivot_wider(names_from = timepoint, values_from = mark) %>%
    column_to_rownames("gene") %>%
    as.matrix()
  
  # Define symmetric color scale around 0
  max_abs = max(abs(mat_heat), na.rm = TRUE)
  breaks = seq(-max_abs, max_abs, length.out = 101)
  
  # Plot heatmap
  p2 = pheatmap(mat_heat,
           cluster_cols = FALSE,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           breaks = breaks,
           display_numbers = overlay_mat,
           number_color = "black",
           fontsize_number = 10,
           show_rownames = TRUE,
           fontsize_row = 10,
           cellwidth = 25,
           cellheight = 12,
           main = "DEGs: * = LFC>2 & pval<=0.05, X = pval>0.05")
  ggsave(paste0("figures/volcano_heatmap_plots/", trt, "_Holm_4fold.png"), p2, width = 13.2, height = 13.2)
}
  