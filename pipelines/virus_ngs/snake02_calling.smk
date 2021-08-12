import pandas as pd
import numpy as np
from snakemake.utils import validate
from snakemake.utils import min_version
from snakemake.shell import shell

# min_version("5.18.0")
# min_version("5.7.1")

# report: "../report/workflow.rst"
# container: "continuumio/miniconda3:4.8.2"


###### Config file and sample sheets #####
# config
configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

# process sample table - QC - selecting samples with genome coverage cutoff
perc_genome = config['params']['coverage']['percent_genome_cov']
at_depth = config['params']['coverage']['at_depth']

samples = pd.read_csv(config["samples_QC"])
samples = samples.query("{0} > {1}".format(at_depth, perc_genome))

samples.set_index("sample", drop=False, inplace=True)

validate(samples, schema="schemas/samples.schema.yaml")

##### Target rules #####
PILON = expand("pilon/{sample}.fasta", sample=samples.index)

# helper functions
include: "rules/common.smk"

rule all:
    input: "variants/combined/all.filt.PASS.GTset.vcf.gz",
        PILON

##### Modules #####
include: "rules/calling.smk"
include: "rules/assembly.smk"
