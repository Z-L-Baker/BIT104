#!/bin/bash

##useful shortcuts
stardir=/mnt/scratch/c1120287/BK_RNA_Seq/star_mapping
genomedir=/mnt/scratch/c1120287/hs_genome
export stardir
export genomedir

trimdir=/mnt/scratch/c1120287/BK_RNA_Seq/trim_reads
export trimdir

for f in ${trimdir}/*_1.fastq.gz 
do
R1=$(basename $f | cut -f1 -d.)
base=$(echo $R1 | sed 's/trim_1/trim/')

export base

sbatch 3-star-loop.sh

done
