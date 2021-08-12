#!/bin/bash
export PATH="/usr/local/bin/miniconda3/bin:$PATH"
source activate /usr/local/bin/miniconda3/envs/qiime2-2020.11
export PYTHONPATH=/usr/local/src:/usr/local/bin/miniconda3/envs/qiime2-2020.11/lib/python3.6/site-packages:/usr/local/lib/python3.7/dist-packages:/usr/lib/python3/dist-packages
script_full_path=$(dirname "$0")
python3 $script_full_path/qiime2_16s.py "$@"
