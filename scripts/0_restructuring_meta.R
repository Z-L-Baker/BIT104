## Script to restructure the meta data in order to combine control groups into one super-control (a choice based on evidence from prior analysis)

#### Set-up ####
rm(list=ls())

workDir = "/mnt/clusters/admiral/data/c21082179/BIT104/New_Analysis"  # working directory for the R session
setwd(workDir)

# Load in the metadata
meta = read.table("metadata/metadata_complete.txt", header=TRUE, sep="\t", na.strings="")
str(meta)

# Combine the carrier and water into the combined control group 
# basically if either Carrier or Water are in the treatment column replace them with Control (ID and filename will keep the speciic info about them for when required)
combined_meta = meta
combined_meta$Treatment = ifelse(meta$Treatment %in% c("Carrier", "Water"), "Control", meta$Treatment)

# Re-create the Treatment time column to reflect this change
combined_meta$Timetreatment = paste0(combined_meta$Treatment, combined_meta$Time)

# Save this combined meta table to be called in scripts further down the pipeline
write.table(combined_meta, "metadata/combined_control_metadata.txt", row.names= F, sep ="\t", quote = F)
