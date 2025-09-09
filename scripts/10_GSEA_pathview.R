#### GSEA and Pathview R SCRIPT ####

#### Set-up ####
rm(list=ls())
#Set working directory
workDir = "/mnt/clusters/admiral/data/c21082179/BIT104/New_Analysis"  # working directory for the R session
setwd(workDir)

library(dplyr)
library(clusterProfiler)
library(pathview)

## Read in the .RData file from the omnibus_method_results.R script :>
load("Filtered_DESeq2_Results.RData")

#make kegg pathway name to ID map
hsa_kegg = download_KEGG("hsa")
path_name_to_id = setNames(hsa_kegg$KEGGPATHID2NAME$from, hsa_kegg$KEGGPATHID2NAME$to)

# load kegg sets with entrez ids (the gene universe)
# and the translation df
kegg_sets = readRDS("sets/kegg_sets_entrez.rds")
kegg_sets_fixed = stack(kegg_sets)
colnames(kegg_sets_fixed) = c("gene", "term")
head(kegg_sets_fixed)
kegg_df = kegg_sets_fixed[, c(2, 1)]

translate = read.table("sets/kegg_gene2entrez.tsv", header=TRUE, sep="\t", na.strings="")
row.names(translate) = translate$geneID

# create results directory
dir.create("figures/GSEA_pathview")

for(resnum in 1:6){
  dir.create(paste0("figures/GSEA_pathview/res", resnum))
  df_name = paste0("res", resnum, "_BH_table")
  # Use results table
  gsea_dat = get(df_name)
  # Genes that didn't pass the screening are not signifcant but still need to be in the universe for gsea so I will set them to 0 (no change)
  gsea_dat = gsea_dat %>% replace(is.na(.), 0)

  # translate genes into entrez ids (this will drop genes without a translation)
  gsea_dat = merge(gsea_dat, translate, by = "row.names")
  gsea_dat = gsea_dat %>% dplyr::select(-Row.names, -geneID)
  row.names(gsea_dat) = gsea_dat$entrezgene_id

  # keep only the LFC columns
  lfc_mat = gsea_dat %>% dplyr::select(dplyr::matches("shrk_L2FC"))  
  colnames(lfc_mat) = sub("_shrk_L2FC$", "", colnames(lfc_mat))
  
  ranked_lists = lapply(colnames(lfc_mat), function(col){
    vec = lfc_mat[[col]]
    names(vec) = rownames(lfc_mat)
  
    # sort decreasing
    sort(vec, decreasing = TRUE)
  })


  names(ranked_lists) = colnames(lfc_mat)

  gsea_results = setNames(
  lapply(ranked_lists, function(geneVec) {
      GSEA(
        geneList      = geneVec,
        TERM2GENE     = kegg_df,
        minGSSize     = 10,
        maxGSSize     = 500,
        # keep all pathways for global FDR correction across hypotheses
        pvalueCutoff  = 1,
        pAdjustMethod = "BH",
        by = "fgsea"
      )
    }),
    names(ranked_lists)
  )

  # extract all pvals for global FDR
  all_pvals = unlist(lapply(gsea_results, function(res) res$pvalue))
  global_fdr = p.adjust(all_pvals, method = "BH")
  counter = 1
  for (i in seq_along(gsea_results)) {
    n = nrow(gsea_results[[i]]@result)
    gsea_results[[i]]@result$padj_global = global_fdr[counter:(counter + n - 1)]
    counter = counter + n
  }
  
  # do the pcut off by the global padj
  gsea_results_sig = lapply(gsea_results, function(res) {
    res@result = res@result[res@result$padj_global < 0.05, ]
    res
  })
  
  for(i in 1:length(gsea_results_sig)){
    res_df = as.data.frame(gsea_results_sig[[i]])
    hypo = names(gsea_results_sig)[[i]]
    
    dir.create(paste0("figures/GSEA_pathview/res", resnum, "/", hypo))
    write.table(res_df, paste0("figures/GSEA_pathview/res", resnum, "/", hypo, "/GSEA_sigres.tsv"), row.names= F, sep ="\t", quote = F)
    
    sig_paths = res_df$ID
  
    #pathview only works for KEGG IDs... not names :< 
  
    # map significant pathways
    sig_path_ids = path_name_to_id[sig_paths]
    
    # because why would pathview NOT have an option for output directory. can't fathom why I might not want to output hundreds of .png files directly into my workDir
    outDir = paste0("figures/GSEA_pathview/res", resnum, "/", hypo)
    setwd(outDir)
    
    for(pid in sig_path_ids){
      
      ### SOME PIDS THROW ERRORS WITH PATHVIEW.... No current fix I can find anywhere... will skip identified issue ones
      if(pid %in% c("hsa01230")){
        message("Skipping pathway: ", pid, " due to issues with pathview formatting...")
        next
      } else {
      pathview(
        gene.data  = ranked_lists[[i]],
        pathway.id = pid,
        species    = "hsa",
        same.layer = FALSE
      )}
    }
    
    setwd(workDir)
  }
}