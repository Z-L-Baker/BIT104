rm(list=ls())

workDir = "/mnt/clusters/admiral/data/c21082179/BIT104/transcriptome_testing"  # working directory for the R session
setwd(workDir)

clstr_files = list.files(path = "nr_transcriptomes/", pattern = "\\.clstr$", full.names = TRUE)

get_cluster_sizes = function(filepath) {
  lines = readLines(filepath)
  sizes = c()
  count = 0
  for (line in lines) {
    if (startsWith(line, ">Cluster")) {
      if (count > 0) sizes = c(sizes, count)
      count = 0
    } else {
      count = count + 1
    }
  }
  if (count > 0) sizes = c(sizes, count)
  return(sizes)
}

library(tidyverse)

cluster_data = map_dfr(clstr_files, function(file) {
  sizes = get_cluster_sizes(file)
  data.frame(
    size = sizes,
    identity = basename(file)
  )
})

head(cluster_data)

cluster_data %>%
  group_by(identity) %>%
  summarise(
    num_clusters = n(),
    mean_size = mean(size),
    median_size = median(size),
    max_size = max(size),
    min_size = min(size),
    singletons = sum(size == 1),
    singleton_percent = round(100 * sum(size == 1) / n(), 2)
  )

ggplot(cluster_data, aes(x = size)) +
  geom_histogram(binwidth = 1, fill = "steelblue") +
  facet_wrap(~ identity) +
  labs(title = "Cluster Size Distributions by Threshold", x = "Cluster Size", y = "Frequency")
