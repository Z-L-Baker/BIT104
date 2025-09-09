#!/bin/bash
#SBATCH --partition=defq      # the requested queue
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

# Add new software module path

export MODULEPATH=/trinity/shared/apps/spack/share/spack/modules/linux-rocky8-ivybridge:$MODULEPATH

# Load subread version compiled on jumbo

module load subread/2.0.6-qefcod7

##define varibles
#star map output dir
mkupdir=/mnt/clusters/admiral/data/c21082179/BIT104/aligning/markup

#genome dir
genomedir=/mnt/clusters/admiral/data/c21082179/BIT104/aligning/genome

#featurecount
featuredir=/mnt/clusters/admiral/data/c21082179/BIT104/aligning/featurecount

#trimdir
trimdir=/mnt/clusters/admiral/data/c21082179/BIT104/aligning/trimdata

#define array

file=$(ls ${trimdir}/*_R1.fastq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)

R1=$(basename $file | cut -f1 -d.)
base=$(echo $R1 | sed 's/_trim_R1$//')

featureCounts -T "${SLURM_CPUS_PER_TASK}" -F GTF -t exon -g gene_id \
	-a "${genomedir}/validated.gtf" \
	-o "${featuredir}/${base}.markdup.featurecount" \
	"${mkupdir}/${base}.markdup.bam"

featureCounts -T ${SLURM_CPUS_PER_TASK} -F GTF -t exon -g gene_id \
	-a "${genomedir}/validated.gtf" \
	-o "${featuredir}/${base}.rmdup.featurecount" \
	"${mkupdir}/${base}.rmdup.bam"

