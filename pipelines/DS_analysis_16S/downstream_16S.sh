#!/bin/bash
export PATH="/usr/local/bin/miniconda3/bin:$PATH"
export PYTHONUNBUFFERED=1
source /usr/local/bin/miniconda3/etc/profile.d/conda.sh
conda activate /home/admin/.conda/envs/qiime2-2018.6
export PYTHONPATH=/usr/local/src:/usr/local/lib/python3.5/dist-packages/
script_full_path=$(dirname "$0")
python3 $script_full_path/downstream_16S.py "$@"
