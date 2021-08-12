import pandas as pd
import numpy as np
from snakemake.utils import validate
from snakemake.utils import min_version
from snakemake.shell import shell

# min_version("5.18.0")
# min_version("5.7.1")

# report: "../report/workflow.rst"
# container: "continuumio/miniconda3:4.8.2"

# function to add proper paths to sample metadata
def add_fq(row):
    read_path = "data/raw/"
    if row['library'] == 'PAIRED':
        return ["{}{}_1.fastq".format(read_path, row['sample']), "{}{}_2.fastq".format(read_path, row['sample'])]
    else:
        return ["{}{}.fastq".format(read_path, row['sample']), np.nan]

###### Config file and sample sheets #####
# config
configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

# process sample table
samples = pd.read_csv(config["samples"])
samples['fq1'], samples['fq2'] = zip(*samples.apply(lambda row: add_fq(row), axis=1))
samples.set_index("sample", drop=False, inplace=True)

validate(samples, schema="schemas/samples.schema.yaml")

##### Target rules #####

# helper functions
include: "rules/common.smk"

rule all:
    input: 
        config["samples_QC"],
        "qc/multiqc_report.html"

##### Modules #####
include: "rules/mapping.smk"
include: "rules/qc.smk"
