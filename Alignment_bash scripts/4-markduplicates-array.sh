#!/bin/bash
#SBATCH --partition=defq      # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     #
#SBATCH --cpus-per-task=2     #
#SBATCH --mem-per-cpu=7000     # in megabytes, unless unit explicitly stated
#SBATCH --error=%J.err         # redirect stderr to this file
#SBATCH --output=%J.out        # redirect stdout to this file
#SBATCH --array=1-324%81
##SBATCH --mail-user=[insert email address]@Cardiff.ac.uk  # email address used for event notification
##SBATCH --mail-type=all                                   # email on job start, failure and end


echo "Some Usable Environment Variables:"
echo "================================="
echo "hostname=$(hostname)"
echo \$SLURM_JOB_ID=${SLURM_JOB_ID}
echo \$SLURM_NTASKS=${SLURM_NTASKS}
echo \$SLURM_NTASKS_PER_NODE=${SLURM_NTASKS_PER_NODE}
echo \$SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
echo \$SLURM_JOB_CPUS_PER_NODE=${SLURM_JOB_CPUS_PER_NODE}
echo \$SLURM_MEM_PER_CPU=${SLURM_MEM_PER_CPU}

# Write jobscript to output file (good for reproducability)
cat $0

##define varibles
#star map output dir
stardir=/mnt/clusters/admiral/data/c21082179/BIT104/aligning/star

#genome dir
mkupdir=/mnt/clusters/admiral/data/c21082179/BIT104/aligning/markup

#trimdir
trimdir=/mnt/clusters/admiral/data/c21082179/BIT104/aligning/trimdata

#define array

file=$(ls ${trimdir}/*_R1.fastq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)

R1=$(basename $file | cut -f1 -d.)
base=$(echo $R1 | sed 's/_trim_R1$//')

#load some modules
module load picard/3.1.1-dyntalv
module load samtools/1.19.2-i77kweo

samtools index ${stardir}/${base}-unsort.Aligned.out.bam

samtools sort -@ ${SLURM_CPUS_PER_TASK} -o ${stardir}/${base}.sorted.bam ${stardir}/${base}-unsort.Aligned.out.bam

##  MARK DUPLICATES  ##
picard MarkDuplicates I=${stardir}/${base}.sorted.bam \
	              O=${mkupdir}/${base}.markdup.bam \
		      M=${mkupdir}/${base}.metrics.markdup.txt \
		      REMOVE_DUPLICATES=false \
		      VALIDATION_STRINGENCY=SILENT

## REMOVE DUPLICATES ##
picard MarkDuplicates I=${stardir}/${base}.sorted.bam \
	              O=${mkupdir}/${base}.rmdup.bam \
		      M=${mkupdir}/${base}.metrics.rmdup.txt \
		      REMOVE_DUPLICATES=true \
		      VALIDATION_STRINGENCY=SILENT

