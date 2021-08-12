#!/bin/bash

# submit this file with:  qsub submit.sh

# run from current working directory
#$ -cwd

#$ -m be
#$ -M brendan.jeffrey@nih.gov

# log dirs
#$ -e ./log/submit_log/
#$ -o ./log/submit_log/

snake_log=$PWD/log/snake_log
mkdir -p $snake_log

# create qsub command
sbcmd="qsub -l {cluster.mem} -e $snake_log -o $snake_log "

module load snakemake || exit 1

# mapping reads, QC
snakemake -R multiqc --snakefile snake01_mapping.smk --use-conda --jobs 100 --rerun-incomplete --keep-going --cluster-config cluster.yaml --cluster "$sbcmd" --latency-wait 120 all

# calling and pilon assembly
snakemake --snakefile snake02_calling.smk --use-conda --jobs 100 --rerun-incomplete --keep-going --cluster-config cluster.yaml --cluster "$sbcmd" --latency-wait 120 all
