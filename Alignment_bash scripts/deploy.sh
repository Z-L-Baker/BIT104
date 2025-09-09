#!/bin/bash

cd /mnt/clusters/admiral/data/c21082179/BIT104/aligning

sbatch --job-name=test_job -d singleton scripts/1A-QC-array.sh
sbatch --job-name=test_job -d singleton scripts/1B-trim-array.sh
sbatch --job-name=test_job -d singleton scripts/1C-QC-array.sh
#sbatch --job-name=test_gtf -d singleton scripts/2-star_index_genome.sh
sbatch --job-name=test_job -d singleton scripts/3-star-array.sh
sbatch --job-name=test_job -d singleton scripts/4-markduplicates-array.sh
sbatch --job-name=test_job -d singleton scripts/5-featurecounts-array.sh
