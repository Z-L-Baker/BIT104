#!/bin/bash
#SBATCH --partition=defq     # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     #
#SBATCH --cpus-per-task=4     #
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

##trim reads

module load fastp/0.23.4-mhegxff

cp "${rawdir}/${base}_R1.fastq.gz" "/tmp/${base}_R1.fastq.gz"
#cp "${rawdir}/${base}_R2.fastq.gz" "/tmp/${base}_R2.fastq.gz"

# Run fastp
fastp -w ${SLURM_CPUS_PER_TASK} \
      -i "/tmp/${base}_R1.fastq.gz" \
      -o "/tmp/${base}_trim_R1.fastq.gz" \
      --trim_tail1 20 \
      --trim_poly_x \
      --length_required 30 \
      --low_complexity_filter \
      --complexity_threshold 30 \
      -h "${trimdir}/posttrim_QC/${base}_fastp.html" \
      -j "${trimdir}/posttrim_QC/${base}_fastp.json"

cp "/tmp/${base}_trim_R1.fastq.gz" "${trimdir}/${base}_trim_R1.fastq.gz"
#cp "/tmp/${base}_trim_R2.fastq.gz" "${trimdir}/${base}_trim_R2.fastq.gz"

rm "/tmp/${base}_"*".fastq.gz"

module unload fastp/0.23.4-mhegxff
