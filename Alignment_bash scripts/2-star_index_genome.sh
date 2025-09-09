#!/bin/bash
#SBATCH --partition=jumbo      # the requested queue
#SBATCH --nodes=1              # number of nodes to use
#SBATCH --tasks-per-node=1     #
#SBATCH --cpus-per-task=32     #
#SBATCH --mem=512000	       # in megabytes, unless unit explicitly stated
#SBATCH --error=%J.err         # redirect stderr to this file
#SBATCH --output=%J.out        # redirect stdout to this file
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

# Load some modules
module load STAR/2.7.6a
 
## Useful shortcuts
export refdir=/mnt/scratch2/GROUP-smbpk/DSS_timecourse/genome

## Change --sjdbOverhang to length of your sequence data /2 minus 1

STAR --runThreadN ${SLURM_JOB_CPUS_PER_NODE} \
        --limitGenomeGenerateRAM 512000000000 \
	--runMode genomeGenerate \
	--genomeDir  $refdir/ \
	--genomeFastaFiles $refdir/GCA_905160935.1_enchytraeus_crypticus_complete_genome_genomic.fna \
	--sjdbGTFfile $refdir/E_crypticus_genome_annotation.gtf \
	--sjdbOverhang 49
