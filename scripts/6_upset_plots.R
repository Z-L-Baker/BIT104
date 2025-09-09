#### R SCRIPT FOR MAKING UPSET PLOTS FROM DESEQ2 OUTPUT DATA ####

#### Set-up ####
rm(list=ls())
#Set working directory
workDir = "/mnt/clusters/admiral/data/c21082179/BIT104/New_Analysis"  # working directory for the R session
setwd(workDir)

#install.packages("UpSetR")

# Load libraries
library(UpSetR)
library("grid")
library("ggplotify")
library("gifski")


dir.create("figures/upsets")
dir.create("figures/upsets/treatments")

#### Treatment Time-Averaged Upset Plots ####
fp = "DEG_lists/TreatmentAverage_BH_2fold_shrk/all_DEGs"
files = list.files(fp, pattern = "\\.txt$", full.names = T)

dir.create("figures/upsets/treatments/averaged")

DEG_lists = lapply(files, readLines)
names(DEG_lists) = basename(files)
names(DEG_lists) = sub("_.*", "", names(DEG_lists))

Azo_Pro_Group = DEG_lists[c("Azoxy", "Azo-Pro", "Pro")]
Cyp_Pro_Group = DEG_lists[c("Cyp", "Cyp-Pro", "Pro")]
Imd_Pro_Group = DEG_lists[c("Imid", "Imd-Pro", "Pro")]

ALL_avg = upset(fromList(DEG_lists),
               keep.order = T,
               order.by = "freq",
               main.bar.color="steelblue", 
               mainbar.y.label = "Time-Averaged Treatment Intersections", 
               sets.x.label = "Number of DEGs", 
               text.scale = 2, 
               empty.intersections = T)
AP_avg = upset(fromList(Azo_Pro_Group),
               sets = c( "Pro", "Azo-Pro", "Azoxy"),
               keep.order = T,
               order.by = "degree",
               main.bar.color="#FF9933", 
               mainbar.y.label = "Time-Averaged Treatment Intersections",
               sets.x.label = "Number of DEGs", 
               text.scale = 2, 
               empty.intersections = T)
CP_avg = upset(fromList(Cyp_Pro_Group), 
               sets = c( "Pro", "Cyp-Pro", "Cyp"),
               keep.order = T,
               order.by = "degree",
               main.bar.color="#00CC33", 
               mainbar.y.label = "Time-Averaged Treatment Intersections",
               sets.x.label = "Number of DEGs", 
               text.scale = 2, 
               empty.intersections = T)
IP_avg = upset(fromList(Imd_Pro_Group),
               sets = c( "Pro", "Imd-Pro", "Imid"),
               keep.order = T,
               order.by = "degree",
               main.bar.color="#33CCFF", 
               mainbar.y.label = "Time-Averaged Treatment Intersections",
               sets.x.label = "Number of DEGs", 
               text.scale = 2, 
               empty.intersections = T)

ALL_avg = as.ggplot(ALL_avg, hjust = -0.1)
AP_avg = as.ggplot(AP_avg, hjust = -0.1)
CP_avg = as.ggplot(CP_avg, hjust = -0.1)
IP_avg = as.ggplot(IP_avg, hjust = -0.1)

ggsave("figures/upsets/treatments/averaged/All_avg.png", ALL_avg)
ggsave("figures/upsets/treatments/averaged/Azo_Pro_avg.png", AP_avg)
ggsave("figures/upsets/treatments/averaged/Cyp_Pro_avg.png", CP_avg)
ggsave("figures/upsets/treatments/averaged/Imd_Pro_avg.png", IP_avg)

#### Treatment at 1hour Upset Plots ####

fp = "DEG_lists/TreatmentPerTime_BH_4fold_shrk/all_DEGs"
files = list.files(fp, pattern = "\\vs_Control_1.txt$", full.names = T)

DEG_lists = lapply(files, readLines)
names(DEG_lists) = basename(files)
names(DEG_lists) = sub("_.*", "", names(DEG_lists))

Azo_Pro_Group = DEG_lists[c("Azoxy", "Azo-Pro", "Pro")]
Cyp_Pro_Group = DEG_lists[c("Cyp", "Cyp-Pro", "Pro")]
Imd_Pro_Group = DEG_lists[c("Imid", "Imd-Pro", "Pro")]

ALL_1 = upset(fromList(DEG_lists),
                keep.order = T,
                order.by = "freq",
                main.bar.color="steelblue", 
                mainbar.y.label = "Treatment Intersections at 1 hour",
                #mainbar.y.max = 90,
                sets.x.label = "Number of DEGs", 
                text.scale = 2, 
                empty.intersections = T)
AP_1 = upset(fromList(Azo_Pro_Group),
             sets = c( "Pro", "Azo-Pro", "Azoxy"),
             keep.order = T,
             order.by = "degree",
             main.bar.color="#FF9933", 
             mainbar.y.label = "Treatment Intersections at 1 hour",
             #mainbar.y.max = 125,
             sets.x.label = "Number of DEGs", 
             text.scale = 2, 
             empty.intersections = T)
CP_1 = upset(fromList(Cyp_Pro_Group),
             sets = c( "Pro", "Cyp-Pro", "Cyp"),
             keep.order = T,
             order.by = "degree",
             main.bar.color="#00CC33", 
             mainbar.y.label = "Treatment Intersections at 1 hour",
             #mainbar.y.max = 95,
             sets.x.label = "Number of DEGs", 
             text.scale = 2, 
             empty.intersections = T)
IP_1 = upset(fromList(Imd_Pro_Group),
             sets = c( "Pro", "Imd-Pro", "Imid"),
             keep.order = T,
             order.by = "degree",
             main.bar.color="#33CCFF", 
             mainbar.y.label = "Treatment Intersections at 1 hour",
             #mainbar.y.max = 115,
             sets.x.label = "Number of DEGs", 
             text.scale = 2, 
             empty.intersections = T)

ALL_1 = as.ggplot(ALL_1, hjust = -0.1)
AP_1 = as.ggplot(AP_1, hjust = -0.1)
CP_1 = as.ggplot(CP_1, hjust = -0.1)
IP_1 = as.ggplot(IP_1, hjust = -0.1)

ggsave("figures/upsets/treatments/All_1.png", ALL_1)
ggsave("figures/upsets/treatments/Azo_Pro_1.png", AP_1)
ggsave("figures/upsets/treatments/Cyp_Pro_1.png", CP_1)
ggsave("figures/upsets/treatments/Imd_Pro_1.png", IP_1)

#### Treatment at 4hour Upset Plots ####

fp = "DEG_lists/TreatmentPerTime_BH_4fold_shrk/all_DEGs"
files = list.files(fp, pattern = "\\vs_Control_4.txt$", full.names = T)

DEG_lists = lapply(files, readLines)
names(DEG_lists) = basename(files)
names(DEG_lists) = sub("_.*", "", names(DEG_lists))

Azo_Pro_Group = DEG_lists[c("Azoxy", "Azo-Pro", "Pro")]
Cyp_Pro_Group = DEG_lists[c("Cyp", "Cyp-Pro", "Pro")]
Imd_Pro_Group = DEG_lists[c("Imid", "Imd-Pro", "Pro")]

ALL_4 = upset(fromList(DEG_lists),
              keep.order = T,
              order.by = "freq",
              main.bar.color="steelblue", 
              mainbar.y.label = "Treatment Intersections at 4 hours",
              #mainbar.y.max = 90,
              sets.x.label = "Number of DEGs", 
              text.scale = 2, 
              empty.intersections = T)
AP_4 = upset(fromList(Azo_Pro_Group),
             sets = c( "Pro", "Azo-Pro", "Azoxy"),
             keep.order = T,
             order.by = "degree",
             main.bar.color="#FF9933", 
             mainbar.y.label = "Treatment Intersections at 4 hours",
             #mainbar.y.max = 125,
             sets.x.label = "Number of DEGs", 
             text.scale = 2, 
             empty.intersections = T)
CP_4 = upset(fromList(Cyp_Pro_Group),
             sets = c( "Pro", "Cyp-Pro", "Cyp"),
             keep.order = T,
             order.by = "degree",
             main.bar.color="#00CC33", 
             mainbar.y.label = "Treatment Intersections at 4 hours",
             #mainbar.y.max = 95,
             sets.x.label = "Number of DEGs", 
             text.scale = 2, 
             empty.intersections = T)
IP_4 = upset(fromList(Imd_Pro_Group),
             sets = c( "Pro", "Imd-Pro", "Imid"),
             keep.order = T,
             order.by = "degree",
             main.bar.color="#33CCFF", 
             mainbar.y.label = "Treatment Intersections at 4 hours",
             #mainbar.y.max = 115,
             sets.x.label = "Number of DEGs", 
             text.scale = 2, 
             empty.intersections = T)

ALL_4 = as.ggplot(ALL_4, hjust = -0.1)
AP_4 = as.ggplot(AP_4, hjust = -0.1)
CP_4 = as.ggplot(CP_4, hjust = -0.1)
IP_4 = as.ggplot(IP_4, hjust = -0.1)

ggsave("figures/upsets/treatments/All_4.png", ALL_4)
ggsave("figures/upsets/treatments/Azo_Pro_4.png", AP_4)
ggsave("figures/upsets/treatments/Cyp_Pro_4.png", CP_4)
ggsave("figures/upsets/treatments/Imd_Pro_4.png", IP_4)


#### Treatment at 8hour Upset Plots ####

fp = "DEG_lists/TreatmentPerTime_BH_4fold_shrk/all_DEGs"
files = list.files(fp, pattern = "\\vs_Control_8.txt$", full.names = T)

DEG_lists = lapply(files, readLines)
names(DEG_lists) = basename(files)
names(DEG_lists) = sub("_.*", "", names(DEG_lists))

Azo_Pro_Group = DEG_lists[c("Azoxy", "Azo-Pro", "Pro")]
Cyp_Pro_Group = DEG_lists[c("Cyp", "Cyp-Pro", "Pro")]
Imd_Pro_Group = DEG_lists[c("Imid", "Imd-Pro", "Pro")]

ALL_8 = upset(fromList(DEG_lists),
              keep.order = T,
              order.by = "freq",
              main.bar.color="steelblue", 
              mainbar.y.label = "Treatment Intersections at 8 hours",
              #mainbar.y.max = 90,
              sets.x.label = "Number of DEGs", 
              text.scale = 2, 
              empty.intersections = T)
AP_8 = upset(fromList(Azo_Pro_Group),
             sets = c( "Pro", "Azo-Pro", "Azoxy"),
             keep.order = T,
             order.by = "degree",
             main.bar.color="#FF9933", 
             mainbar.y.label = "Treatment Intersections at 8 hours",
             #mainbar.y.max = 125,
             sets.x.label = "Number of DEGs", 
             text.scale = 2, 
             empty.intersections = T)
CP_8 = upset(fromList(Cyp_Pro_Group),
             sets = c( "Pro", "Cyp-Pro", "Cyp"),
             keep.order = T,
             order.by = "degree",
             main.bar.color="#00CC33", 
             mainbar.y.label = "Treatment Intersections at 8 hours",
             #mainbar.y.max = 95,
             sets.x.label = "Number of DEGs", 
             text.scale = 2, 
             empty.intersections = T)
IP_8 = upset(fromList(Imd_Pro_Group),
             sets = c( "Pro", "Imd-Pro", "Imid"),
             keep.order = T,
             order.by = "degree",
             main.bar.color="#33CCFF", 
             mainbar.y.label = "Treatment Intersections at 8 hours",
             #mainbar.y.max = 115,
             sets.x.label = "Number of DEGs", 
             text.scale = 2, 
             empty.intersections = T)

ALL_8 = as.ggplot(ALL_8, hjust = -0.1)
AP_8 = as.ggplot(AP_8, hjust = -0.1)
CP_8 = as.ggplot(CP_8, hjust = -0.1)
IP_8 = as.ggplot(IP_8, hjust = -0.1)

ggsave("figures/upsets/treatments/All_8.png", ALL_8)
ggsave("figures/upsets/treatments/Azo_Pro_8.png", AP_8)
ggsave("figures/upsets/treatments/Cyp_Pro_8.png", CP_8)
ggsave("figures/upsets/treatments/Imd_Pro_8.png", IP_8)


#### Treatment at 12hour Upset Plots ####

fp = "DEG_lists/TreatmentPerTime_BH_4fold_shrk/all_DEGs"
files = list.files(fp, pattern = "\\vs_Control_12.txt$", full.names = T)

DEG_lists = lapply(files, readLines)
names(DEG_lists) = basename(files)
names(DEG_lists) = sub("_.*", "", names(DEG_lists))

Azo_Pro_Group = DEG_lists[c("Azoxy", "Azo-Pro", "Pro")]
Cyp_Pro_Group = DEG_lists[c("Cyp", "Cyp-Pro", "Pro")]
Imd_Pro_Group = DEG_lists[c("Imid", "Imd-Pro", "Pro")]

ALL_12 = upset(fromList(DEG_lists),
              keep.order = T,
              order.by = "freq",
              main.bar.color="steelblue", 
              mainbar.y.label = "Treatment Intersections at 12 hours",
              #mainbar.y.max = 90,
              sets.x.label = "Number of DEGs", 
              text.scale = 2, 
              empty.intersections = T)
AP_12 = upset(fromList(Azo_Pro_Group),
             sets = c( "Pro", "Azo-Pro", "Azoxy"),
             keep.order = T,
             order.by = "degree",
             main.bar.color="#FF9933", 
             mainbar.y.label = "Treatment Intersections at 12 hours",
             #mainbar.y.max = 125,
             sets.x.label = "Number of DEGs", 
             text.scale = 2, 
             empty.intersections = T)
CP_12 = upset(fromList(Cyp_Pro_Group),
             sets = c( "Pro", "Cyp-Pro", "Cyp"),
             keep.order = T,
             order.by = "degree",
             main.bar.color="#00CC33", 
             mainbar.y.label = "Treatment Intersections at 12 hours",
             #mainbar.y.max = 95,
             sets.x.label = "Number of DEGs", 
             text.scale = 2, 
             empty.intersections = T)
IP_12 = upset(fromList(Imd_Pro_Group),
             sets = c( "Pro", "Imd-Pro", "Imid"),
             keep.order = T,
             order.by = "degree",
             main.bar.color="#33CCFF", 
             mainbar.y.label = "Treatment Intersections at 12 hours",
             #mainbar.y.max = 115,
             sets.x.label = "Number of DEGs", 
             text.scale = 2, 
             empty.intersections = T)

ALL_12 = as.ggplot(ALL_12, hjust = -0.1)
AP_12 = as.ggplot(AP_12, hjust = -0.1)
CP_12 = as.ggplot(CP_12, hjust = -0.1)
IP_12 = as.ggplot(IP_12, hjust = -0.1)

ggsave("figures/upsets/treatments/All_12.png", ALL_12)
ggsave("figures/upsets/treatments/Azo_Pro_12.png", AP_12)
ggsave("figures/upsets/treatments/Cyp_Pro_12.png", CP_12)
ggsave("figures/upsets/treatments/Imd_Pro_12.png", IP_12)


#### Treatment at 16hour Upset Plots ####

fp = "DEG_lists/TreatmentPerTime_BH_4fold_shrk/all_DEGs"
files = list.files(fp, pattern = "\\vs_Control_16.txt$", full.names = T)

DEG_lists = lapply(files, readLines)
names(DEG_lists) = basename(files)
names(DEG_lists) = sub("_.*", "", names(DEG_lists))

Azo_Pro_Group = DEG_lists[c("Azoxy", "Azo-Pro", "Pro")]
Cyp_Pro_Group = DEG_lists[c("Cyp", "Cyp-Pro", "Pro")]
Imd_Pro_Group = DEG_lists[c("Imid", "Imd-Pro", "Pro")]

ALL_16 = upset(fromList(DEG_lists),
               keep.order = T,
               order.by = "freq",
               main.bar.color="steelblue", 
               mainbar.y.label = "Treatment Intersections at 16 hours",
               #mainbar.y.max = 90,
               sets.x.label = "Number of DEGs", 
               text.scale = 2, 
               empty.intersections = T)
AP_16 = upset(fromList(Azo_Pro_Group),
              sets = c( "Pro", "Azo-Pro", "Azoxy"),
              keep.order = T,
              order.by = "degree",
              main.bar.color="#FF9933", 
              mainbar.y.label = "Treatment Intersections at 16 hours",
              #mainbar.y.max = 125,
              sets.x.label = "Number of DEGs", 
              text.scale = 2, 
              empty.intersections = T)
CP_16 = upset(fromList(Cyp_Pro_Group),
              sets = c( "Pro", "Cyp-Pro", "Cyp"),
              keep.order = T,
              order.by = "degree",
              main.bar.color="#00CC33", 
              mainbar.y.label = "Treatment Intersections at 16 hours",
              #mainbar.y.max = 95,
              sets.x.label = "Number of DEGs", 
              text.scale = 2, 
              empty.intersections = T)
IP_16 = upset(fromList(Imd_Pro_Group),
              sets = c( "Pro", "Imd-Pro", "Imid"),
              keep.order = T,
              order.by = "degree",
              main.bar.color="#33CCFF", 
              mainbar.y.label = "Treatment Intersections at 16 hours",
              #mainbar.y.max = 115,
              sets.x.label = "Number of DEGs", 
              text.scale = 2, 
              empty.intersections = T)

ALL_16 = as.ggplot(ALL_16, hjust = -0.1)
AP_16 = as.ggplot(AP_16, hjust = -0.1)
CP_16 = as.ggplot(CP_16, hjust = -0.1)
IP_16 = as.ggplot(IP_16, hjust = -0.1)

ggsave("figures/upsets/treatments/All_16.png", ALL_16)
ggsave("figures/upsets/treatments/Azo_Pro_16.png", AP_16)
ggsave("figures/upsets/treatments/Cyp_Pro_16.png", CP_16)
ggsave("figures/upsets/treatments/Imd_Pro_16.png", IP_16)


#### Treatment at 24hour Upset Plots ####

fp = "DEG_lists/TreatmentPerTime_BH_4fold_shrk/all_DEGs"
files = list.files(fp, pattern = "\\vs_Control_24.txt$", full.names = T)

DEG_lists = lapply(files, readLines)
names(DEG_lists) = basename(files)
names(DEG_lists) = sub("_.*", "", names(DEG_lists))

Azo_Pro_Group = DEG_lists[c("Azoxy", "Azo-Pro", "Pro")]
Cyp_Pro_Group = DEG_lists[c("Cyp", "Cyp-Pro", "Pro")]
Imd_Pro_Group = DEG_lists[c("Imid", "Imd-Pro", "Pro")]

ALL_24 = upset(fromList(DEG_lists),
               keep.order = T,
               order.by = "freq",
               main.bar.color="steelblue", 
               mainbar.y.label = "Treatment Intersections at 24 hours",
               #mainbar.y.max = 90,
               sets.x.label = "Number of DEGs", 
               text.scale = 2, 
               empty.intersections = T)
AP_24 = upset(fromList(Azo_Pro_Group),
              sets = c( "Pro", "Azo-Pro", "Azoxy"),
              keep.order = T,
              order.by = "degree",
              main.bar.color="#FF9933", 
              mainbar.y.label = "Treatment Intersections at 24 hours",
              #mainbar.y.max = 125,
              sets.x.label = "Number of DEGs", 
              text.scale = 2, 
              empty.intersections = T)
CP_24 = upset(fromList(Cyp_Pro_Group),
              sets = c( "Pro", "Cyp-Pro", "Cyp"),
              keep.order = T,
              order.by = "degree",
              main.bar.color="#00CC33", 
              mainbar.y.label = "Treatment Intersections at 24 hours",
              #mainbar.y.max = 95,
              sets.x.label = "Number of DEGs", 
              text.scale = 2, 
              empty.intersections = T)
IP_24 = upset(fromList(Imd_Pro_Group),
              sets = c( "Pro", "Imd-Pro", "Imid"),
              keep.order = T,
              order.by = "degree",
              main.bar.color="#33CCFF", 
              mainbar.y.label = "Treatment Intersections at 24 hours",
              #mainbar.y.max = 115,
              sets.x.label = "Number of DEGs", 
              text.scale = 2, 
              empty.intersections = T)

ALL_24 = as.ggplot(ALL_24, hjust = -0.1)
AP_24 = as.ggplot(AP_24, hjust = -0.1)
CP_24 = as.ggplot(CP_24, hjust = -0.1)
IP_24 = as.ggplot(IP_24, hjust = -0.1)

ggsave("figures/upsets/treatments/All_24.png", ALL_24)
ggsave("figures/upsets/treatments/Azo_Pro_24.png", AP_24)
ggsave("figures/upsets/treatments/Cyp_Pro_24.png", CP_24)
ggsave("figures/upsets/treatments/Imd_Pro_24.png", IP_24)


#### ALL Group TC GIFs ####
#png_files = list.files(
#  "figures/upsets/treatments/",
#  pattern = "^All.*\\.png$",
#  full.names = TRUE
#)
#nums = as.numeric(gsub(".*All_([0-9]+)\\.png$", "\\1", basename(png_files)))
#png_files = png_files[order(nums)]
#
#gifski(png_files, gif_file = "figures/upsets/treatments/All_upset.gif", delay = 2)  # delay in seconds per frame
#
##### AP Group TC GIFs ####
#png_files = list.files(
#  "figures/upsets/treatments/",
#  pattern = "^Azo_Pro.*\\.png$",
#  full.names = TRUE
#)
#nums = as.numeric(gsub(".*Azo_Pro_([0-9]+)\\.png$", "\\1", basename(png_files)))
#png_files = png_files[order(nums)]
#
#gifski(png_files, gif_file = "figures/upsets/treatments/Azo_Pro_upset.gif", delay = 2)  # delay in seconds per frame
#
##### CP Group TC GIFs ####
#png_files = list.files(
#  "figures/upsets/treatments/",
#  pattern = "^Cyp_Pro.*\\.png$",
#  full.names = TRUE
#)
#nums = as.numeric(gsub(".*Cyp_Pro_([0-9]+)\\.png$", "\\1", basename(png_files)))
#png_files = png_files[order(nums)]
#
#gifski(png_files, gif_file = "figures/upsets/treatments/Cyp_Pro_upset.gif", delay = 2)  # delay in seconds per frame
#
##### IP Group TC GIFs ####
#png_files = list.files(
#  "figures/upsets/treatments/",
#  pattern = "^Imd_Pro.*\\.png$",
#  full.names = TRUE
#)
#nums = as.numeric(gsub(".*Imd_Pro_([0-9]+)\\.png$", "\\1", basename(png_files)))
#png_files = png_files[order(nums)]
#
#gifski(png_files, gif_file = "figures/upsets/treatments/Imd_Pro_upset.gif", delay = 2)  # delay in seconds per frame


#### Intreraction Time-Averaged Upset Plots ####
dir.create("figures/upsets/interactions")
fp = "DEG_lists/InteractionAverage_BH_2fold_shrk/all_DEGs"
files = list.files(fp, pattern = "\\.txt$", full.names = T)


DEG_lists = lapply(files, readLines)
names(DEG_lists) = basename(files)
names(DEG_lists) = sub(".txt$", "", names(DEG_lists))


Avg = upset(fromList(DEG_lists),
                keep.order = T,
                order.by = "freq",
                main.bar.color="steelblue", 
                mainbar.y.label = "Time-Averaged Interaction Intersections", 
                sets.x.label = "Number of DEGs", 
                text.scale = 2, 
                empty.intersections = T)


Avg = as.ggplot(Avg, hjust = -0.1)


ggsave("figures/upsets/interactions/Avg.png", Avg)

#### Interaction at 1hour Upset Plots ####

fp = "DEG_lists/InteractionPerTime_BH_4fold_shrk/all_DEGs"
files = list.files(fp, pattern = "\\_1h.txt$", full.names = T)

DEG_lists = lapply(files, readLines)
names(DEG_lists) = basename(files)
names(DEG_lists) = sub(".txt$", "", names(DEG_lists))

Int_1 = upset(fromList(DEG_lists),
              keep.order = T,
              order.by = "freq",
              main.bar.color="steelblue", 
              mainbar.y.label = "Interaction Intersections at 1 hour",
              #mainbar.y.max = 70,
              sets.x.label = "Number of DEGs", 
              text.scale = 2, 
              empty.intersections = T)


Int_1 = as.ggplot(Int_1, hjust = -0.1)

ggsave("figures/upsets/interactions/Int_1.png", Int_1)


#### Interaction at 4hour Upset Plots ####

fp = "DEG_lists/InteractionPerTime_BH_4fold_shrk/all_DEGs"
files = list.files(fp, pattern = "\\_4h.txt$", full.names = T)

DEG_lists = lapply(files, readLines)
names(DEG_lists) = basename(files)
names(DEG_lists) = sub(".txt$", "", names(DEG_lists))

Int_4 = upset(fromList(DEG_lists),
              keep.order = T,
              order.by = "freq",
              main.bar.color="steelblue", 
              mainbar.y.label = "Interaction Intersections at 4 hours",
              #mainbar.y.max = 70,
              sets.x.label = "Number of DEGs", 
              text.scale = 2, 
              empty.intersections = T)


Int_4 = as.ggplot(Int_4, hjust = -0.1)

ggsave("figures/upsets/interactions/Int_4.png", Int_4)

#### Interaction at 8hour Upset Plots ####

fp = "DEG_lists/InteractionPerTime_BH_4fold_shrk/all_DEGs"
files = list.files(fp, pattern = "\\_8h.txt$", full.names = T)

DEG_lists = lapply(files, readLines)
names(DEG_lists) = basename(files)
names(DEG_lists) = sub(".txt$", "", names(DEG_lists))

Int_8 = upset(fromList(DEG_lists),
              keep.order = T,
              order.by = "freq",
              main.bar.color="steelblue", 
              mainbar.y.label = "Interaction Intersections at 8 hours",
              #mainbar.y.max = 70,
              sets.x.label = "Number of DEGs", 
              text.scale = 2, 
              empty.intersections = T)


Int_8 = as.ggplot(Int_8, hjust = -0.1)

ggsave("figures/upsets/interactions/Int_8.png", Int_8)

#### Interaction at 12hour Upset Plots ####

fp = "DEG_lists/InteractionPerTime_BH_4fold_shrk/all_DEGs"
files = list.files(fp, pattern = "\\_12h.txt$", full.names = T)

DEG_lists = lapply(files, readLines)
names(DEG_lists) = basename(files)
names(DEG_lists) = sub(".txt$", "", names(DEG_lists))

Int_12 = upset(fromList(DEG_lists),
              keep.order = T,
              order.by = "freq",
              main.bar.color="steelblue", 
              mainbar.y.label = "Interaction Intersections at 12 hours",
              #mainbar.y.max = 70,
              sets.x.label = "Number of DEGs", 
              text.scale = 2, 
              empty.intersections = T)


Int_12 = as.ggplot(Int_12, hjust = -0.1)

ggsave("figures/upsets/interactions/Int_12.png", Int_12)


#### Interaction at 16hour Upset Plots ####

fp = "DEG_lists/InteractionPerTime_BH_4fold_shrk/all_DEGs"
files = list.files(fp, pattern = "\\_16h.txt$", full.names = T)

DEG_lists = lapply(files, readLines)
names(DEG_lists) = basename(files)
names(DEG_lists) = sub(".txt$", "", names(DEG_lists))

Int_16 = upset(fromList(DEG_lists),
               keep.order = T,
               order.by = "freq",
               main.bar.color="steelblue", 
               mainbar.y.label = "Interaction Intersections at 16 hours",
               #mainbar.y.max = 70,
               sets.x.label = "Number of DEGs", 
               text.scale = 2, 
               empty.intersections = T)


Int_16 = as.ggplot(Int_16, hjust = -0.1)

ggsave("figures/upsets/interactions/Int_16.png", Int_16)

#### Interaction at 24hour Upset Plots ####

fp = "DEG_lists/InteractionPerTime_BH_4fold_shrk/all_DEGs"
files = list.files(fp, pattern = "\\_24h.txt$", full.names = T)

DEG_lists = lapply(files, readLines)
names(DEG_lists) = basename(files)
names(DEG_lists) = sub(".txt$", "", names(DEG_lists))

Int_24 = upset(fromList(DEG_lists),
               keep.order = T,
               order.by = "freq",
               main.bar.color="steelblue", 
               mainbar.y.label = "Interaction Intersections at 24 hours",
               #mainbar.y.max = 70,
               sets.x.label = "Number of DEGs", 
               text.scale = 2, 
               empty.intersections = T)


Int_24 = as.ggplot(Int_24, hjust = -0.1)

ggsave("figures/upsets/interactions/Int_24.png", Int_24)

#### ALL Group TC GIFs ####
#png_files = list.files(
#  "figures/upsets/interactions/",
#  pattern = "^Int.*\\.png$",
#  full.names = TRUE
#)
#nums = as.numeric(gsub(".*Int_([0-9]+)\\.png$", "\\1", basename(png_files)))
#png_files = png_files[order(nums)]

#gifski(png_files, gif_file = "figures/upsets/interactions/Int_upset.gif", delay = 2)  # delay in seconds per frame


