#### COMPARITOR OF KEGG AND GO SETS ####

#### Set-up ####
rm(list=ls())
#Set working directory
workDir = "/mnt/clusters/admiral/data/c21082179/BIT104/New_Analysis"  # working directory for the R session
setwd(workDir)

kegg_sets = readRDS("sets/kegg_sets_final.rds")
go_sets = readRDS("sets/go_sets.rds")

both_sets_genes = intersect(
  unique(unlist(go_sets)),
  unique(unlist(kegg_sets))
)

length(unique(unlist(go_sets)))
length(unique(unlist(kegg_sets)))
length(unique(unlist(both_sets_genes)))
# 3628 genes appear in both sets :< bu regardless, this will allow us to 

# there were 9439 possible unique genes across the actual human kegg sets...
# 3682 of these genes were found among the final 7556 unique genes with corresponding entrez IDs from the data... and thus in the kegg_sets 
# (this data will have collapsed paralogs)
# there were 12989 unique genes which had corresponding GO terms and were thus sorted into go term clusters
# (this data should have all paralogs)

# Being able to back up KEGG findings with this data would be immensely helpful because obviously the KEGG results are based on such a small subset
# if 3649 genes appear in both sets that means almost all of the KEGG relevant genes found in the data appear in the GO sets...
# upon filtering down to these (so as not to bias results due the comparably HUGE go sets) could find overlaps
# thus support KEGG pathways and bring more clarity to vague go clusters if they align well.

go_sets_filt   = lapply(go_sets,   function(gs) intersect(gs, both_sets_genes))
kegg_sets_filt = lapply(kegg_sets, function(ks) intersect(ks, both_sets_genes))

to_binary_matrix = function(gene_sets, allgenes) {
  mat = sapply(gene_sets, function(genes) allgenes %in% genes)
  rownames(mat) = allgenes
  return(mat)
}

# Made into binary matrices to be faster to process
go_mat = to_binary_matrix(go_sets_filt, both_sets_genes)
kegg_mat = to_binary_matrix(kegg_sets_filt, both_sets_genes)

# results list for all possible combinations of kegg pathways and go clusters
results_list = vector("list", length = ncol(go_mat) * ncol(kegg_mat))
# counter to track which combination is being tested
counter = 1
# if there aren't at LEAST 5 genes overlapping between groups throw the whole combo out to save time
min_overlap = 5

# run loop, for each go cluter compare each kegg cluster and then move onto the next
for(i in seq_len(ncol(go_mat))) {
  go_vec = go_mat[, i]
  
  for(j in seq_len(ncol(kegg_mat))) {
    kegg_vec = kegg_mat[, j]
    
    # find the genes in both vectors
    intersection = sum(go_vec & kegg_vec)
    # find the genes in either vector
    union = sum(go_vec | kegg_vec)
    # calculate Jaccard's index
    jaccard = intersection / union
    
    # skip pairs below minimum overlap
    if(intersection < min_overlap) next
    
    # make a 2x2 grid for Fisher's exact test
    # genes in both
    a = intersection
    # ones in go but not kegg
    b = sum(go_vec & !kegg_vec)
    # ones in kegg but not go
    c = sum(!go_vec & kegg_vec)
    # neither
    d = sum(!go_vec & !kegg_vec)
    
    #calculare fishers and then pull the p val and the odds ratio
    fisher_table = matrix(c(a, b, c, d), nrow = 2)
    fisher_res = fisher.test(fisher_table, alternative = "greater")
    fisher_p = fisher_res$p.value
    fisher_odds = fisher_res$estimate
    
    # save the results to the results list
    results_list[[counter]] = data.frame(
      GO_set = colnames(go_mat)[i],
      KEGG_set = colnames(kegg_mat)[j],
      overlap = intersection,
      jaccard = jaccard,
      fisher_odds = fisher_odds,
      fisher_p = fisher_p
    )
    # increase counter (aka move onto next combo)
    counter = counter + 1
  }
}

# make results a df
results_df = bind_rows(results_list)

# FDR mitigation
results_df$fisher_p_adj = p.adjust(results_df$fisher_p, method = "BH")

# now filter for actual pairings that are any good.
# HAS to be significant in fishers, with strong odds ratio
# HAS also got to be BOTH actually got some overlap AND have that overlap be like 20% of the overall set size.
sig_results = results_df %>%
  filter(
    fisher_p_adj < 0.05,
    overlap >= 5,
    jaccard >= 0.2,
    fisher_odds >= 4
  )

#save results
write.table(sig_results, "sets/significant_overlaps.tsv", row.names= F, sep ="\t", quote = F)
