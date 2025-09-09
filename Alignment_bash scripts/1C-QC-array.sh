#!/bin/bash
#SBATCH --partition=defq     # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     #
#SBATCH --cpus-per-task=1     #
#SBATCH --mem-per-cpu=1000     # in megabytes, unless unit explicitly stated
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
#rawdata
rawdir=/mnt/clusters/admiral/data/c21082179/BIT104/aligning/rawdata

#trimdir
trimdir=/mnt/clusters/admiral/data/c21082179/BIT104/aligning/trimdata

#define array

file=$(ls ${rawdir}/*_R1.fastq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)

R1=$(basename $file | cut -f1 -d.)
base=$(echo $R1 | sed 's/_R1$//')


##Fastqc trimmed reads

module load fastqc/v0.12.1

# fastqc the outputs
fastqc -t ${SLURM_CPUS_PER_TASK} "${trimdir}/${base}_trim_R1.fastq.gz" --outdir "${trimdir}/posttrim_QC"
#fastqc -t ${SLURM_CPUS_PER_TASK} "${trimdir}/${base}_trim_R2.fastq.gz" --outdir "${trimdir}/posttrim_QC"


module unload fastqc/v0.12.1
