rm(list=ls())

workDir = "/mnt/clusters/admiral/data/c21082179/BIT104"  # working directory for the R session
setwd(workDir)

rsem = read.table("transcript_abundance/Azo_Pro_12hour_1/RSEM.genes.results", header=TRUE)
fc = read.table("Analysis/data/featurecount_20bp/Azo_Pro_12hour_1.markdup.featurecount", sep="\t", quote="\"", header=T, skip=0, stringsAsFactors=FALSE)
rsem = rsem[,c(1, 5)]
fc = fc[,c(1, 7)]
colnames(rsem) = c("GeneID", "RSEM")
colnames(fc) = c("GeneID", "FeatureCounts")

hist(log2(rsem$RSEM + 1), main="RSEM", xlab="log2(count + 1)")
hist(log2(fc$FeatureCounts + 1), main="FeatureCounts", xlab="log2(count + 1)")

sum(rsem$RSEM)
sum(fc$FeatureCounts)

sum(rsem$RSEM > 0)
sum(fc$FeatureCounts > 0)

sum(rsem$RSEM == 0)
sum(fc$FeatureCounts == 0)

par(mfrow=c(1,2))
hist(log2(rsem$RSEM + 1), main="RSEM", xlab="log2(count + 1)", col="skyblue", breaks=50)
hist(log2(fc$FeatureCounts + 1), main="FeatureCounts", xlab="log2(count + 1)", col="salmon", breaks=50)

library(tximport)
txi.rsem = tximport("transcript_abundance/Azo_Pro_1hour_12/RSEM.genes.results", type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem$counts)

hist(log2(txi.rsem$counts + 1), main="RSEM", xlab="log2(count + 1)")
hist(log2(fc$FeatureCounts + 1), main="FeatureCounts", xlab="log2(count + 1)")

sum(txi.rsem$counts)
sum(fc$FeatureCounts)

sum(txi.rsem$counts > 0)
sum(fc$FeatureCounts > 0)

sum(txi.rsem$counts == 0)
sum(fc$FeatureCounts == 0)

par(mfrow=c(1,2))
hist(log2(txi.rsem$counts + 1), main="RSEM", xlab="log2(count + 1)", col="skyblue", breaks=50)
hist(log2(fc$FeatureCounts + 1), main="FeatureCounts", xlab="log2(count + 1)", col="salmon", breaks=50)




log_rsem = log2(txi.rsem$counts + 1)
log_fc   = log2(fc$FeatureCounts + 1)

log_rsem = log2(txi.rsem$counts[txi.rsem$counts > 0] + 1)
log_fc   = log2(fc$FeatureCounts[fc$FeatureCounts > 0] + 1)


xlim_range = range(c(log_rsem, log_fc))
ylim_range = range( 
  0, 
  hist(log_rsem, plot = FALSE, breaks = 50)$counts,
  hist(log_fc, plot = FALSE, breaks = 50)$counts
)

# Plot side by side
par(mfrow = c(1, 2))  # 1 row, 2 columns

hist(log_rsem,
     main = "RSEM",
     xlab = "log2(count + 1)",
     col = "skyblue",
     breaks = 50,
     xlim = xlim_range,
     ylim = ylim_range
)

hist(log_fc,
     main = "FeatureCounts",
     xlab = "log2(count + 1)",
     col = "salmon",
     breaks = 50,
     xlim = xlim_range,
     ylim = ylim_range
)

par(mfrow = c(1, 1))  # reset layout
