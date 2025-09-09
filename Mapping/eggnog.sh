#!/bin/bash
#SBATCH --partition=longq     
#SBATCH --nodes=1              
#SBATCH --tasks-per-node=1     
#SBATCH --cpus-per-task=16      
#SBATCH --mem=32G     

echo "Some Usable Environment Variables:"
echo "================================="
echo "hostname=$(hostname)"
echo \$SLURM_JOB_ID=${SLURM_JOB_ID}
echo \$SLURM_NTASKS=${SLURM_NTASKS}
echo \$SLURM_NTASKS_PER_NODE=${SLURM_NTASKS_PER_NODE}
echo \$SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
echo \$SLURM_JOB_CPUS_PER_NODE=${SLURM_JOB_CPUS_PER_NODE}
echo \$SLURM_MEM_PER_CPU=${SLURM_MEM_PER_CPU}

cat $0

cd /mnt/clusters/admiral/data/c21082179/BIT104/eggnog

pipedir=$(pwd)

# load singularity module
module load apptainer

SINGULARITY_IMAGE_NAME=eggnog-mapper.2.1.13.sif

if [ -f ${pipedir}/singularities/${SINGULARITY_IMAGE_NAME} ]; then
    echo "Singularity image exists"
else
    echo "Singularity image does not exist"
        singularity pull eggnog-mapper.2.1.13.sif \
        https://depot.galaxyproject.org/singularity/eggnog-mapper%3A2.1.13--pyhdfd78af_1
fi

echo ${singularities}

SINGIMAGEDIR=${pipedir}/singularities
SINGIMAGENAME=${SINGULARITY_IMAGE_NAME}

# Set working directory
WORKDIR=${pipedir}
SING_IMAGE=${WORKDIR}/singularities/eggnog-mapper.2.1.13.sif

export BINDS="${BINDS},${WORKDIR}:${WORKDIR}"

# Sanity check: list emapper.py inside container
apptainer exec --bind ${WORKDIR}:${WORKDIR} ${SING_IMAGE} ls -l /usr/local/bin/emapper.py

apptainer exec --bind ${WORKDIR}:${WORKDIR} ${SING_IMAGE} python3 --version

#apptainer exec --bind ${WORKDIR}:${WORKDIR} ${SING_IMAGE} download_eggnog_data.py -F -y --data_dir ${pipedir}/eggnog_data

## Run emapper.py
apptainer exec singularities/eggnog-mapper.2.1.13.sif \
	emapper.py \
	-i ${WORKDIR}/E_crypticus_protein.fasta \
        -o ${WORKDIR}/ec_eggnog \
	--data_dir ${WORKDIR}/eggnog_data \
	--override \
	> ${WORKDIR}/ec.full.log 2> ${WORKDIR}/ec.err.log
