#!/bin/bash

# Job Name
#$ -N test
#$ -cwd
#$ -m be
#$ -M your.email@nih.gov
#$ -j y
#$ -o log/

# create directory that will contain split SRA ids file if doesnt exist
temp_sras='../temp_sra_ids'
mkdir -p $temp_sras

# main sample table
sra_table='../samples_dev10.csv'
sra_ids='../sra_ids.txt'

# generate just ids file
awk -F"," 'NR>1 {print $1}' $sra_table > $sra_ids

# need number of lines in each subfile
lines=10
split --numeric-suffixes=1 --lines=$lines --suffix-length=3 $sra_ids $temp_sras/sra.1