#! /bin/bash
#SBATCH -o /home/rapleeid/practice/nextflow/sysout/nextflow%j.txt
#SBATCH -e /home/rapleeid/practice/nextflow/sysout/nextflow-er%j.txt

set -e

. /data/rapleeid/conda/etc/profile.d/conda.sh
conda activate

module load nextflow
nextflow run main.nf -with-dag /home/rapleeid/practice/nextflow/flowchart.html
#nextflow -resume main.nf
#nextflow help
