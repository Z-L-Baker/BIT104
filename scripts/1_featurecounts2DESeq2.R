## Script to read featurecounts output into DESeq2 - partially adapted from SARtools

# This script takes a few hours to run with my data so it's advisable to submit the script as a slurm job or run in the background with a nohup command - don't try to run it actively and wait.

#### Set-up ####
rm(list=ls())

#BiocManager::install("DESeq2", dependencies=TRUE, lib="/mnt/clusters/admiral/data/c21082179/containers/Rpackages")
.libPaths("/mnt/clusters/admiral/data/c21082179/containers/Rpackages")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DESeq2")

library(DESeq2)


workDir = "/mnt/clusters/admiral/data/c21082179/BIT104/New_Analysis"  # working directory for the R session
setwd(workDir)

#### Load in the metadata ####
meta = read.table("metadata/combined_control_metadata.txt", header=TRUE, sep="\t", na.strings="")

str(meta)

# Set the variants of interest as factors
meta$Treatment = as.factor(meta$Treatment)
meta$Time = as.factor(meta$Time)
meta$Timetreatment = as.factor(meta$Timetreatment)


#### Importing count data ####
# ADAPTED FROM SARTOOLS LOADCOUNTDATA FUNCTION: https://github.com/PF2-pasteur-fr/SARTools/blob/master/R/loadCountData.R
# Only a few changes made to explicitly adjust for my data
  
# Assign metadata info to varibles to search for IDs and associated files
labels = as.character(meta[,1])
files = as.character(meta[,2])

# Search for the first file name within the given directory and read that file in
rawCounts = read.table(file.path("data/featurecount", files[1]), sep="\t", quote="\"", header=T, skip=0, stringsAsFactors=FALSE)
# This is featurecounts output so the gene ID is in column 1 and the counts in column 7 and this is all the data we need for DESeq2 so drop everything else
rawCounts = rawCounts[,c(1, 7)]
#Rename the first column to be Id with the corresponding ID for that file as outlined in the metadata
colnames(rawCounts) = c("Id", labels[1])
# Check for any duplicate gene IDs in the file and stop if there are any - they will cause problems downstream
if (any(duplicated(rawCounts$Id))){
  stop("Duplicated feature names in ", files[1], ": ", 
        paste(unique(rawCounts$Id[duplicated(rawCounts$Id)]), collapse=", "))
}
# Print a summary of the file, how many genes/features are present and how many of these have 0 counts
cat("Loading files:\n")
cat(files[1],": ",length(rawCounts[,labels[1]])," rows and ",sum(rawCounts[,labels[1]]==0)," null count(s)\n",sep="")
  
# Do this same process for ALL the files iteratively with a temp dataframe and then merge them into the rawCounts dataframe
for (i in 2:length(files)){
  tmp = read.table(file.path("data/featurecount", files[i]), sep="\t", quote="\"", header=T, skip=0, stringsAsFactors=FALSE)
  tmp = tmp[,c(1, 7)]
  colnames(tmp) = c("Id", labels[i])
  if (any(duplicated(tmp$Id))){
    stop("Duplicated feature names in ", files[i], ": ", 
          paste(unique(tmp$Id[duplicated(tmp$Id)]), collapse=", "))
  }
  rawCounts = merge(rawCounts, tmp, by="Id", all=TRUE)
  cat(files[i],": ",length(tmp[,labels[i]])," rows and ",sum(tmp[,labels[i]]==0)," null count(s)\n",sep="")
}

# Replace any NA values with 0
rawCounts[is.na(rawCounts)] = 0
# DESeq2 requires the count data as a matrix, so we will drop the gene ID column and convert the rest of the dataframe to a matrix
counts = as.matrix(rawCounts[,-1])
# Then we can assign the gene IDs as the rownames of the new matrix
rownames(counts) = rawCounts[,1]
# Order gene IDs by name
counts = counts[order(rownames(counts)),]

# check that input counts are integers to fit edgeR and DESeq2 requirements
if (any(counts %% 1 != 0)) stop("Input counts are not integer values as required by DESeq2")
  
# Check that the results are reasonable
cat("\nTop of the counts matrix:\n")
print(head(counts))
cat("\nBottom of the counts matrix:\n")
print(tail(counts))
  
  
# Remove any genes where there are ZERO counts across ALL samples because they're unhelpful for scaling methods downstream
counts = counts[rowSums(counts) > 0, ]

#### Flag outliers function ####
# One gene with insane counts in ONE sample where the other biological replicates are all between 0 and 5 counts keeps failing to converge
# function will drop any genes which for their replicate group, have a sample over 100 times!! the median counts of that replicate group. Should protect normal variation.
flag_outlier_genes = function(cts, group) {
  outlier = logical(nrow(cts))
  # iterate through each gene
  for (i in 1:nrow(cts)) {
    # iterate through each replicate group
    for (g in unique(group)) {
      # colelct group counts
      grp_counts = cts[i, group == g]
      # calculate median
      med = median(grp_counts)
      # if any counts in the group are more than 100 times the median
      if (any(grp_counts > 100 * max(1, med))) {
        # label it an outlier :<
        outlier[i] = TRUE
        break
      }
    }
  }
  # return outlier flags
  return(outlier)
}

#### DESeq2 Time-course Analysis Method (Treatment as a factor) ####
  
# This was largely adapted from https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
# Edited and annotated though, as well as adapted for an omnibus testing method.
  
# Importing the clean data as a matrix, adding the model design with the variables treatment and time, and the interaction term
ddsTC = DESeqDataSetFromMatrix(counts, meta, ~ Treatment + Time + Treatment:Time)
  
# Check levels
levels(ddsTC$Treatment)
# Azo-Pro is currently the default level (because it is first alphabetically) so I'm going to change that because I want to compare the treatments against the carrier.
ddsTC$Treatment = relevel(ddsTC$Treatment, ref = "Control") # Note I initially ran this with "Carrier" as the control first, but found minimal differences compared with "Water" so have combined into "Control"
levels(ddsTC$Treatment)  

# 1 hour is first in time though so we don't need to change anything here
levels(ddsTC$Time)

# going to make sure that when there are counts, there are at least 3 samples of 10 
# This SHOULD keep reasonably expressed genes but ensures at least half of a set of replicates actually show the expression 
# (allowing for biological variation between replicates)
ddsTC = estimateSizeFactors(ddsTC)
nc = counts(ddsTC, normalized = TRUE)
keep = rowSums(nc >= 10) >= 3
ddsTC = ddsTC[keep, ]

# One gene with ridiculously high counts in one sample is failing to converge... 
# group samples into biologcial replicates
group = factor(paste(meta$Treatment, meta$Time, sep="_"))

# use outlier function
outlier_flags = flag_outlier_genes(counts(ddsTC), group)
row.names(ddsTC[outlier_flags, ])
ddsTC = ddsTC[!outlier_flags, ]


# Run a likelihood ratio test to determine the significant DEGs of ALL treatment associations (This is one of our omnibus tests to drill down from)
LRT1 = DESeq(ddsTC, test="LRT", reduced = ~ Time)

# Assign the results to a variable to call
resLRT1 = results(LRT1)
  
# Get the gene names from the metadata and put them in the corresponding results for easier identification
resLRT1$symbol = mcols(LRT1)$symbol

# Run ANOTHER likelihood ratio test to determine the significant DEGs of ONLY the treatment time interaction term, another omnibus test to drill down from.
LRT2 = DESeq(ddsTC, test="LRT", reduced = ~ Treatment + Time)

# Assign the results to a variable to call
resLRT2 = results(LRT2)
  
# Get the gene names from the metadata and put them in the corresponding results for easier identification
resLRT2$symbol = mcols(LRT2)$symbol

# Run a wald test we can pull the contrasts from  later - this will be the results we look at for finer details following the omnibus identifications.
WLD = DESeq(ddsTC)

# Assign the results to a variable to call (with only the converged genes)
resWLD = results(WLD)

# Get the gene names from the metadata and put them in the corresponding results for easier identification
resWLD$symbol = mcols(WLD)$symbol



#### DESeq2 Chemical Interactions Analysis Method (Chemical Presence Binaries) ####

# Ok so this will be largely the same but I want to direclty model interaction of treatments so need different chemical presence indicators for Pro, Azo, Cyp and Imid

# Build presence indicators
meta$Azo = factor(ifelse(grepl("^Azo", meta$Treatment), "pres", "abs"))
meta$Cyp = factor(ifelse(grepl("^Cyp", meta$Treatment), "pres", "abs"))
meta$Imid = factor(ifelse(grepl("^Im", meta$Treatment), "pres", "abs"))
meta$Pro = factor(ifelse(grepl("Pro$", meta$Treatment), "pres", "abs"))


# Importing the clean data as a matrix, adding the model design with the variables treatment and time, and the interaction term
ddsInteractions = DESeqDataSetFromMatrix(counts, meta, 
                                         ~ Azo + Cyp + Imid + Pro + Time # Chems
                                         + Azo:Time + Cyp:Time + Imid:Time + Pro:Time # Individual Chem and time combo
                                         + Azo:Pro + Cyp:Pro + Imid:Pro # Chem interactions
                                         + Azo:Pro:Time + Cyp:Pro:Time + Imid:Pro:Time # Chem interactions and time combo
                                         )

# these are all levelled correctly :>

## Need to up max iterations so need to do this manually (normally done through DESeq2 but that doesn't let you increase max iterations for convergence)
ddsInteractions = estimateSizeFactors(ddsInteractions)

# going to make sure that when there are counts, there are at least 3 samples of 10 
# This SHOULD keep reasonably expressed genes but ensures at least half of a set of replicates actually show the expression 
# (allowing for biological variation between replicates)
nc = counts(ddsInteractions, normalized = TRUE)
keep = rowSums(nc >= 10) >= 3
ddsInteractions = ddsInteractions[keep, ]

# One gene with ridiculously high counts in one sample is failing to converge... 
# use outlier function
int_outlier_flags = flag_outlier_genes(counts(ddsInteractions), group)
row.names(ddsInteractions[int_outlier_flags, ])
ddsInteractions = ddsInteractions[!int_outlier_flags, ]

ddsInteractions = estimateDispersions(ddsInteractions)


# Run a likelihood ratio test to determine the significant DEGs of ANY interaction associations (This is one of our omnibus tests to drill down from)
INTLRT1 = nbinomLRT(ddsInteractions, 
                reduced = ~ Azo + Cyp + Imid + Pro + Time # Chems
                + Azo:Time + Cyp:Time + Imid:Time + Pro:Time,# Individual Chem and time combo
                maxit = 500) 

# Assign the results to a variable to call
resINTLRT1 = results(INTLRT1)

# Get the gene names from the metadata and put them in the corresponding results for easier identification
resINTLRT1$symbol = mcols(INTLRT1)$symbol


# Run ANOTHER likelihood ratio test to determine the significant DEGs of ONLY the Chem interactions and time combos, another omnibus test to drill down from.
INTLRT2 = nbinomLRT(ddsInteractions, 
                reduced = ~ Azo + Cyp + Imid + Pro + Time # Chems
                + Azo:Time + Cyp:Time + Imid:Time + Pro:Time # Individual Chem and time combo
                + Azo:Pro + Cyp:Pro + Imid:Pro, # Chem interactions
                maxit = 500)

# Assign the results to a variable to call
resINTLRT2 = results(INTLRT2)

# Get the gene names from the metadata and put them in the corresponding results for easier identification
resINTLRT2$symbol = mcols(INTLRT2)$symbol

# Run a wald test we can pull the contrasts from  later - this will be the results we look at for finer details following the omnibus identifications.
INTWLD = nbinomWaldTest(ddsInteractions, maxit = 500)

# Assign the results to a variable to call (with only the converged genes)
resINTWLD = results(INTWLD)

# Get the gene names from the metadata and put them in the corresponding results for easier identification
resINTWLD$symbol = mcols(INTWLD)$symbol


# Save .RData file
save.image(file=paste0("DESeq2_TC_Results", ".RData"))
