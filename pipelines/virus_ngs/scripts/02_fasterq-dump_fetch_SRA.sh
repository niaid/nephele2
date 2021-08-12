#!/bin/sh

#$ -N fetch_sra
#$ -cwd
#$ -V
#$ -j y
#$ -m e
#$ -M your.email@nih.gov

# usage: qsub -l mem_free=24G 02_fasterq-dump_fetch_SRA.sh

out_reads='../data/raw'
mkdir -p $out_reads
#$ -o $out_reads
#$ -t 1001-1100

# location of sra_ids files
temp_sras='../temp_sra_ids'

# pull reads
module load sra-toolkit

while read SRA; do
    fasterq-dump $SRA --outdir $out_reads
done < $temp_sras/sra.$SGE_TASK_ID

