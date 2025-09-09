#### R SCRIPT FOR MAKING UPSET PLOTS FROM DESEQ2 OUTPUT DATA ####

#### Set-up ####
rm(list=ls())
#Set working directory
workDir = "/mnt/clusters/admiral/data/c21082179/BIT104/New_Analysis"  # working directory for the R session
setwd(workDir)

#install.packages("ComplexUpset")

library(ComplexUpset)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)


dir.create("figures/upsets")
dir.create("figures/upsets/weird/")

fp = "DEG_lists/TreatmentPerTime_BH_4fold_shrk/all_DEGs"
files1 = list.files(fp, pattern = "\\vs_Control_1.txt$", full.names = T)
DEG_lists1 = lapply(files1, readLines)
names(DEG_lists1) = basename(files1)
names(DEG_lists1) = sub("_.*", "", names(DEG_lists1))
files2 = list.files(fp, pattern = "\\vs_Control_4.txt$", full.names = T)
DEG_lists2 = lapply(files2, readLines)
names(DEG_lists2) = basename(files2)
names(DEG_lists2) = sub("_.*", "", names(DEG_lists2))
files3 = list.files(fp, pattern = "\\vs_Control_8.txt$", full.names = T)
DEG_lists3 = lapply(files3, readLines)
names(DEG_lists3) = basename(files3)
names(DEG_lists3) = sub("_.*", "", names(DEG_lists3))
files4 = list.files(fp, pattern = "\\vs_Control_12.txt$", full.names = T)
DEG_lists4 = lapply(files4, readLines)
names(DEG_lists4) = basename(files4)
names(DEG_lists4) = sub("_.*", "", names(DEG_lists4))
files5 = list.files(fp, pattern = "\\vs_Control_16.txt$", full.names = T)
DEG_lists5 = lapply(files5, readLines)
names(DEG_lists5) = basename(files5)
names(DEG_lists5) = sub("_.*", "", names(DEG_lists5))
files6 = list.files(fp, pattern = "\\vs_Control_24.txt$", full.names = T)
DEG_lists6 = lapply(files6, readLines)
names(DEG_lists6) = basename(files6)
names(DEG_lists6) = sub("_.*", "", names(DEG_lists6))

all_timepoints = list(
  t1 = DEG_lists1,
  t4 = DEG_lists2,
  t8 = DEG_lists3,
  t12 = DEG_lists4,
  t16 = DEG_lists5,
  t24 = DEG_lists6
)

df = map_dfr(names(all_timepoints), function(tp) {
  sets = all_timepoints[[tp]]
  map_dfr(names(sets), function(setname) {
    tibble(
      gene = sets[[setname]],
      set  = setname,
      timepoint = tp
    )
  })
})

df_wide = df %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from = set,
    values_from = present,
    values_fill = list(present = 0)
  )

df_azo_set = df_wide %>% select(gene, timepoint, `Azo-Pro`, Azoxy, Pro)
df_cyp_set = df_wide %>% select(gene, timepoint, `Cyp-Pro`, Cyp, Pro)
df_imid_set = df_wide %>% select(gene, timepoint, `Imd-Pro`, Imid, Pro)

sets = list(
  Azo_Pro = df_azo_set,
  Cyp_Pro = df_cyp_set,
  Imid_Pro = df_imid_set)

for(set in names(sets)){
  df = data.frame(sets[[set]])   # use [[ ]] not [ ]
  colnames(df) = gsub(paste0("^", set, "\\."), "", colnames(df))
  colnames(df) = gsub(".Pro$", "_Pro", colnames(df))
  
  names = colnames(df)[3:5]   # set columns
  names2 = colnames(df)[2:5]  # timepoint + sets
  
  intersection_counts = df %>%
    group_by(across(all_of(c("timepoint", names)))) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(intersection = apply(pick(all_of(names)), 1, function(row) {
      sets_on = names(row)[row == 1]
      paste(sets_on, collapse = " ∩ ")
    })) %>%
    filter(intersection != "")
  
  # make sure timepoints are ordered correctly
  intersection_counts = intersection_counts %>%
    mutate(timepoint = factor(timepoint, levels = c("t1","t4","t8","t12","t16","t24")))
  
  plot = ggplot(intersection_counts, aes(x = timepoint, y = count, group = intersection, color = intersection)) +
    geom_line(linewidth = 2) +
    geom_point(size = 1) +
    labs(x = "Time (hours)", y = "DEGs in Intersection", color = "Intersection") +
    theme_minimal(base_size = 20) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(paste0("figures/upsets/weird/",set, "BH_4fold.png"), plot)
}
