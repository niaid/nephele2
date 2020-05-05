#!/bin/bash

## export PATH=/usr/local/bin/miniconda3/bin:$PATH:/usr/local/sbin:/usr/sbin:/sbin  ## for non-root user

export PATH="/usr/local/bin/miniconda3/bin:$PATH"
source activate /home/admin/.conda/envs/qiime2-2018.6
export PYTHONPATH=/usr/local/src:/usr/local/lib/python3.5/dist-packages/
script_full_path=$(dirname "$0")
python3 $script_full_path/qc_reads.py "$@"
