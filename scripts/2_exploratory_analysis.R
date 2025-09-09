## Script to visualise experiment variability AND to perform statistical power analysis for later

rm(list=ls())
#Set working directory
workDir = "/mnt/clusters/admiral/data/c21082179/BIT104/New_Analysis"  # working directory for the R session
setwd(workDir)

## Read in the .RData file from the featurecounts2DESeq2.R script
load("DESeq2_TC_Results.RData")

##Load libraries
library("dplyr")
library("ggplot2")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("tidyr")
library("gganimate")
library("gifski")
library("grid")
library("ggplotify")
library("RNASeqPower")

#### pheatmap save function: ####
# Source: https://gist.github.com/mathzero/a2070a24a6b418740c44a5c023f5c01e #
save_pheatmap = function(x, filename, width=12, height=12){
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  if(grepl(".png",filename)){
    png(filename, width=width, height=height, units = "in", res=300)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  else if(grepl(".pdf",filename)){
    pdf(filename, width=width, height=height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  else{
    print("Filename did not contain '.png' or '.pdf'")
  }
}

#### Analysis ####

## Dispersion plots:

p = plotDispEsts(LRT1)
p = as.ggplot(p)
ggsave(plot = p, filename = "figures/plotDispEsts.png", width = 6, height = 3)

# I want to visualise the variability of the normalised counts so I'm going to run a PCA analysis
## Need to transform the results because otherwise the highest genecounts will have the largest absolute differences between the samples

# Variance Stabilising Transformation (VST):
vsd = vst(ddsTC, blind = TRUE)
head(assay(vsd), 3)
colData(vsd)

# Cannot do rlog in this case because we have more than 30 samples which will make rlog VERY computationally inefficient
# Regularized-Logarithm Transformation (RLOG):
#rld = rlog(ddsTC, blind = TRUE)
#head(assay(rld), 3)


# Plotting methods to take a look
dds = estimateSizeFactors(ddsTC)

df = bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))#,
  #as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] = c("x", "y")  

lvls = c("log2(x + 1)", "vst")#, "rlog")
df$transformation = factor(df$transformation, levels=lvls)

tf_comp = ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation) 

ggsave(plot = tf_comp, filename = "figures/tf_comp.png", width = 6, height = 3)

# VST looks good (trend along the y=x with limited variability)
# I will continue using vst for sample variation investigations
pcaData = plotPCA(vsd, intgroup = c("Treatment", "Time"), returnData = TRUE)
pcaData

percentVar = round(100 * attr(pcaData, "percentVar"))

# Define treatment colors (adjust as needed)
treatment_colors = c(
  "Control" = "black", "Azo-Pro" = "#AE4D66", "Azoxy" = "#FF9933",
  "Cyp" = "#00CC33", "Cyp-Pro" = "#2E6666", "Imd-Pro" = "#4866CC",
  "Imid" = "#33CCFF", "Pro" = "#5C0099"
)

sample_pca = ggplot(pcaData, aes(x = PC1, y = PC2, color = Treatment, shape = Time)) +
  geom_point(size =2) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(values = treatment_colors) +
  coord_fixed() +
  ggtitle("PCA with vst data")

ggsave(plot = sample_pca, filename = "figures/sample_pca.png", width = 6, height = 5)

# Plotted with vst data

#### Sample distance heatmap ####
sampleDists = dist(t(assay(vsd)))
sampleDists


sampleDistMatrix = as.matrix(sampleDists)
rownames(sampleDistMatrix) = paste( vsd$Treatment, vsd$Time, vsd$Rep, sep = " - " )
colnames(sampleDistMatrix) = NULL
colors = colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)
sample_dist = pheatmap(sampleDistMatrix,
                        clustering_distance_rows = sampleDists,
                        clustering_distance_cols = sampleDists,
                        col = colors)

save_pheatmap(sample_dist, filename = "figures/sample_distance.png", width = 48, height = 36)

vsd_split = split(colnames(vsd), vsd$Treatment)

for (i in 1:length(vsd_split)){
  name = names(vsd_split[i])
  vsd_obj = vsd[, vsd_split[[i]]]
  sampleDists = dist(t(assay(vsd_obj)))
  sampleDistMatrix = as.matrix(sampleDists)
  rownames(sampleDistMatrix) = paste( vsd_obj$Treatment, vsd_obj$Time, vsd_obj$Rep, sep = " - " )
  colnames(sampleDistMatrix) = NULL
  colors = colorRampPalette( rev(brewer.pal(9, "Reds")) )(255)
  sample_dist = pheatmap(sampleDistMatrix,
                         clustering_distance_rows = sampleDists,
                         clustering_distance_cols = sampleDists,
                         col = colors)
  
  save_pheatmap(sample_dist, filename = paste0("figures/sample_distance_", name, ".png"), width = 12, height = 9)
}

#### Animated PCA to highlight timecourse changes ####

str(pcaData)

pcaData$Time = as.numeric(levels(pcaData$Time))[pcaData$Time]

pcaData = pcaData %>% arrange(Time)

pcaData$Time = as.factor(pcaData$Time)

mean_points = pcaData %>%
  group_by(group) %>%
  summarise(
    mean_PC1 = mean(PC1),
    mean_PC2 = mean(PC2)
  )

mean_points$new = mean_points$group

mean_points = mean_points %>%
  separate(new,
           into = c("Treatment", "Time"),
           sep = ":")


segments = pcaData %>%
  left_join(mean_points, by = "group") %>%
  mutate(
    x = PC1,
    y = PC2,
    xend = mean_PC1,
    yend = mean_PC2
  )

segments$new = segments$group
segments = segments %>%
  separate(new,
           into = c("Treatment", "Time"),
           sep = ":")

# Define treatment colors (adjust as needed)
treatment_colors = c(
  "Control" = "black", "Azo-Pro" = "#AE4D66", "Azoxy" = "#FF9933",
  "Cyp" = "#00CC33", "Cyp-Pro" = "#2E6666", "Imd-Pro" = "#4866CC",
  "Imid" = "#33CCFF", "Pro" = "#5C0099"
)

p = ggplot() +
  geom_point(data = pcaData, aes(x = PC1, y = PC2, color = Treatment), size = 1.5, alpha = 0.3) +
  geom_point(data = mean_points, aes(x = mean_PC1, y = mean_PC2, color = Treatment), size = 5) +
  geom_segment(data = segments, aes(x = xend, y = yend, xend = x, yend = y, color = Treatment), linewidth = 1.5, alpha = 0.3) +
  scale_color_manual(values = treatment_colors) +
  theme_minimal(base_size = 15) +
  labs(title = "PC1 vs PC2: Replicates and TimeTreatment Means",
       x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"))
p

p_anim = p + transition_states(Time) +
  labs(title = "Time (hours): {closest_state}")

anim_obj = animate(
  p_anim,
  renderer = gifski_renderer(loop = TRUE),
  width = 900, height = 600
)

anim_save("figures/PCA_anim.gif", animation = anim_obj)


### rush job traject pca because I can't put gifs in :<

pcaData$Time = as.numeric(levels(pcaData$Time))[pcaData$Time]
pcaData = pcaData %>% arrange(Time)
pcaData$Time = as.factor(pcaData$Time)

mean_points = pcaData %>%
  group_by(group) %>%
  summarise(
    mean_PC1 = mean(PC1),
    mean_PC2 = mean(PC2)
  )


mean_points$new = mean_points$group
mean_points = mean_points %>%
  separate(new, into = c("Treatment", "Time"), sep = ":") %>%
  mutate(Time = as.numeric(Time)) %>%
  arrange(Treatment, Time)  # order by treatment & time


treatment_colors = c(
  "Control" = "black", "Azo-Pro" = "#AE4D66", "Azoxy" = "#FF9933",
  "Cyp" = "#00CC33", "Cyp-Pro" = "#2E6666", "Imd-Pro" = "#4866CC",
  "Imid" = "#33CCFF", "Pro" = "#5C0099"
)


p_static = ggplot() +
  geom_point(data = pcaData, aes(x = PC1, y = PC2, color = Treatment), size = 1.5, alpha = 0.3) +
  geom_point(data = mean_points, aes(x = mean_PC1, y = mean_PC2, color = Treatment), size = 0.5) +
  geom_path(
    data = mean_points,
    aes(x = mean_PC1, y = mean_PC2, color = Treatment, group = Treatment),  # group is key!
    linewidth = 1.2,
    arrow = arrow(length = unit(0.25, "cm"), type = "closed")  # add arrows
  ) +
  scale_color_manual(values = treatment_colors) +
  theme_minimal(base_size = 15) +
  labs(
    title = "PC1 vs PC2: Treatment Trajectories Over Time",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  )

ggsave("figures/PCA_trajectory_version.png", p_static)



#### Calculating sequencing depth and coefficient of variance for power analysis :> ####
# Here I will collect all the info required for RNASeqPower to see the statistical power we have for detecting DE between any 2 given groups in the dataset
# see: https://bioconductor.org/packages/release/bioc/vignettes/RNASeqPower/inst/doc/samplesize.pdf

# Because we are looking at the hypothetical 2 given groups out of the selection I'm going to take the mean sequencing depth per sample across all groups
# I will also for this same reason be finding the median of the cv within groups
# This will give us a generalised overview of the statistical power we have to detect DE between any 2 given groups in the dataset.
# The results will determine the minimum L2FC that can be detected with a power of 0.8 and alpha of 0.05.

# Mean sequencing depth per sampoe
counts_matrix = counts(ddsTC)
depth_per_sample = colSums(counts_matrix)
mean_depth = mean(depth_per_sample)

# make CV function
mat = assay(vsd)
cv = function(x) sd(x) / mean(x)

# assign groups as treatment and save the number of groups
group = colData(ddsTC)$Treatment

# Have to find minimum now they're different due to control merging
group_sizes = table(group)
n_min = min(group_sizes)
t_num = n_min
#t_num = (ncol(ddsTC)/nlevels(group))

# find CV within treatment groups
cv_by_group = sapply(levels(group), function(g) {
  apply(mat[, group == g], 1, cv)
})

cv_all = as.vector(cv_by_group)  # make a vector so as to take median
cv_mediant = median(cv_all, na.rm = TRUE)

min_effect_treatment_tests = rnapower(depth = (mean_depth/1000000), cv = cv_mediant, alpha = 0.05, power = 0.8, n=t_num)

# cv for timetreatment tests
group = colData(ddsTC)$Timetreatment

# Have to find minimum now they're different due to control merging
group_sizes = table(group)
n_min = min(group_sizes)
tt_num = n_min
#tt_num = (ncol(ddsTC)/nlevels(group))

# CV within groups
cv_by_group = sapply(levels(group), function(g) {
  apply(mat[, group == g], 1, cv)
})

cv_all = as.vector(cv_by_group)
cv_mediantt = median(cv_all, na.rm = TRUE)

min_effect_timetreatment_tests = rnapower(depth = (mean_depth/1000000), cv = cv_mediantt, alpha = 0.05, power = 0.8, n=tt_num)

report = c(
  "RNA-seq Power Analysis Summary",
  "------------------------------",
  paste("Mean sequencing depth per sample:", round(mean_depth, 3)),
  paste("Median CV between treatment groups:", round(cv_mediant, 3)),
  paste("Median CV between treatmenttime groups:", round(cv_mediantt, 3)),
  "",
  paste("Calculated treatment differential analyses effect size (", t_num, "samples per treatment type):", round(min_effect_treatment_tests, 3)),
  paste("Hence, minimum Log2FC which can be detected:", round(log2(min_effect_treatment_tests), 3)),
  "",
  paste("Calculated timetreatment differential analyses effect size (", tt_num, "samples per time within treatment type):", round(min_effect_timetreatment_tests, 3)),
  paste("Hence, minimum Log2FC which can be detected:", round(log2(min_effect_timetreatment_tests), 3)),
  "",
  "Notes:",
  "- CVs computed from VST-transformed counts",
  "- Calculated based on alpha of 0.05 and power of 0.8",
  "- Determined with R Package: RNASeqPower_1.46.0"
)

writeLines(report, "rna_seq_power_analysis_report.txt")


## Ok so importantly, we have found results that imply the AVERAGE we can detect is a certain threshold (4), which will be the threshold i will use
# however for detail and transpoarency i't's important to look at ALL comparisons because that's an average and likely some of them will be underpowered (and some overpowered)
# but overpowered aside we need to specify underpowered results because although the 0.8 power is kind of arbitrary, it's important to be clear WHEN it doesn't apply/meet that power threshold
# basically just warn when things are likely being missed - although our positives will be statistically valid we may get a few more false negatives with some comparisons
cv = function(x) sd(x) / mean(x)

# mean depth per sample
counts_matrix = counts(ddsTC)
depth_per_sample = colSums(counts_matrix)

mat = assay(vsd)

# make treatment groups
treatments = colData(ddsTC)$Treatment
treatment_levels = levels(treatments)

# make function to pull the relevant info for each treatment group
treatment_stats = lapply(treatment_levels, function(g) {
  samples = which(treatments == g)
  list(
    group = g,
    n = length(samples),
    depth = mean(depth_per_sample[samples]),
    cv = median(apply(mat[, samples, drop=FALSE], 1, cv), na.rm=TRUE)
  )
})
treatment_stats = do.call(rbind, lapply(treatment_stats, as.data.frame))

# make a function to perform treatment comparisons on pairs of treatments
treat_pairs = combn(treatment_stats$group, 2, simplify=FALSE)
treatment_results = lapply(treat_pairs, function(pair) {
  g1 = treatment_stats[treatment_stats$group == pair[1], ]
  g2 = treatment_stats[treatment_stats$group == pair[2], ]
  
  depth_mean = mean(c(g1$depth, g2$depth))  # sequencing depth averaged still
  effect = rnapower(
    depth = depth_mean/1000000,
    n = g1$n, n2 = g2$n,
    cv = g1$cv, cv2 = g2$cv,
    alpha = 0.05, power = 0.8
  )
  
  data.frame(group1=g1$group, group2=g2$group, min_effect=effect)
})
treatment_results = do.call(rbind, treatment_results)


cv = function(x) sd(x) / mean(x)

# make time treatment groups
ttreatments  = droplevels(colData(ddsTC)$Timetreatment)
tt_levels = levels(ttreatments)

# make function to pull the stats for each time treatment group
ttreatment_stats = do.call(
  rbind,
  lapply(tt_levels, function(g) {
    samples = which(ttreatments == g)
    if (length(samples) == 0) return(NULL)
    
    trt = unique(as.character(treatments[samples]))
    if (length(trt) != 1) {
      warning(sprintf("Group '%s' maps to multiple treatments: %s", g, paste(trt, collapse=", ")))
      trt = trt[1]
    }
    
    data.frame(
      group     = g,
      treatment = trt,
      n         = length(samples),
      depth     = mean(depth_per_sample[samples]),
      cv        = median(apply(mat[, samples, drop = FALSE], 1, cv), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)

## split the stats dataframe into treatment groups (because I don't want to compare timepoints across treatments)
by_trt = split(ttreatment_stats, ttreatment_stats$treatment)

# make a function to take the dataframe into subsections by the treatment type and then make pairs of timepoints to compare WITHIN those subgroups
res_list = lapply(names(by_trt), function(trt) {
  subdf = by_trt[[trt]]
  
  pairs = combn(subdf$group, 2, simplify = FALSE)
  
  # do all the comparisons and save the statistics, minimum effect detected at 80% power and the power at 4-fold effect to report later
  do.call(rbind, lapply(pairs, function(pair) {
    g1 = subdf[subdf$group == pair[1], ]
    g2 = subdf[subdf$group == pair[2], ]
    
    depth_pair = mean(c(g1$depth, g2$depth))
    
    data.frame(
      treatment  = trt,
      group1     = g1$group,
      group2     = g2$group,
      n1         = g1$n,
      n2         = g2$n,
      cv1        = g1$cv,
      cv2        = g2$cv,
      depth_mean = depth_pair,
      min_effect = rnapower(
        depth = depth_pair / 1e6,
        n     = g1$n,
        n2    = g2$n,
        cv    = g1$cv,
        cv2   = g2$cv,
        alpha = 0.05,
        power = 0.8
      ),
      power_at_4fold = rnapower(
        depth = depth_pair / 1e6,
        n     = g1$n,
        n2    = g2$n,
        cv    = g1$cv,
        cv2   = g2$cv,
        alpha = 0.05,
        effect = 4
      ),
      stringsAsFactors = FALSE
    )
  }))
})

ttreatment_results = do.call(rbind, res_list)


control_df = by_trt[["Control"]]

res_list_vs_control = lapply(setdiff(names(by_trt), "Control"), function(trt) {
  subdf = by_trt[[trt]]
  
  do.call(rbind, lapply(seq_len(nrow(subdf)), function(i) {
    g_trt = subdf[i, ]
    
    # strip treatment prefix to get just the timepoint (digits at the end)
    timepoint = sub("\\D+", "", g_trt$group)
    
    # build expected control group name, e.g. "Control1"
    ctrl_group = control_df[control_df$group == paste0("Control", timepoint), ]
    
    if (nrow(ctrl_group) == 0) {
      warning(sprintf("No control group found for %s (timepoint %s)", g_trt$group, timepoint))
      return(NULL)
    }
    
    depth_pair = mean(c(g_trt$depth, ctrl_group$depth))
    
    data.frame(
      treatment   = trt,
      timepoint   = timepoint,
      group_trt   = g_trt$group,
      group_ctrl  = ctrl_group$group,
      n_trt       = g_trt$n,
      n_ctrl      = ctrl_group$n,
      cv_trt      = g_trt$cv,
      cv_ctrl     = ctrl_group$cv,
      depth_mean  = depth_pair,
      min_effect  = rnapower(
        depth = depth_pair / 1e6,
        n     = g_trt$n,
        n2    = ctrl_group$n,
        cv    = g_trt$cv,
        cv2   = ctrl_group$cv,
        alpha = 0.05,
        power = 0.8
      ),
      power_at_4fold = rnapower(
        depth = depth_pair / 1e6,
        n     = g_trt$n,
        n2    = ctrl_group$n,
        cv    = g_trt$cv,
        cv2   = ctrl_group$cv,
        alpha = 0.05,
        effect = 4
      ),
      stringsAsFactors = FALSE
    )
  }))
})

ttreatment_results_vs_control = do.call(rbind, res_list_vs_control)

# flag anything underpowered for transparency purposes
treatment_results$underpowered = treatment_results$min_effect > 2
ttreatment_results$underpowered = ttreatment_results$min_effect > 4
ttreatment_results_vs_control$underpowered = ttreatment_results_vs_control$min_effect > 4

write.table(treatment_results, "treatment_power_analysis.tsv", row.names= F, sep ="\t", quote = F)
write.table(ttreatment_results, "time_treatment_power_analysis.tsv", row.names= F, sep ="\t", quote = F)
write.table(ttreatment_results_vs_control, "treatmentpertime_power_analysis.tsv", row.names= F, sep ="\t", quote = F)
