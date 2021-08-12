# snakefile to merge VCF files
import sys

shell.prefix("set -eo pipefail; ")

# gather sample names  NOTE:  config file must be specified in the snakemake_submit.sh shell script
SAMPLES = config["samples"]


if SAMPLES is None:
    sys.exit("Sample names are needed")

# build sample name for merged VCF
# merged_sample_name = config["experiment_id"] + str(len(SAMPLES))
merged_sample_name = 'Cneoformans13'

# vt package functions being used
PROCESS_LEVEL = ["DbN", "DbND"]

# set paths
RDS_DIR = config["reads"]
BWA_DIR = config["bwa_aln"]
VCF_DIR = config["vcf"]
VCF_MERGE_DIR = config["vcf_merge"]

# define targets, VCF merge
VCF_MERGE = VCF_MERGE_DIR + merged_sample_name + ".merge.vcf.gz"

# define targets, variants, normalized, decomposed
VARIANTS = VCF_MERGE_DIR + merged_sample_name + ".merge.variants.vcf"
NORM = VCF_MERGE_DIR + merged_sample_name  + ".merge.variants.N.vcf"
MULTIFASTA = VCF_MERGE_DIR + merged_sample_name + ".merge.snps.fa"

# define rules
rule all:
    # input: VCF_MERGE, NORM, MULTIFASTA
    input: NORM, MULTIFASTA

rule generate_merge_list:
    """
    need to generate a file with samples to merge, argument list too long error workaround
    """
    input: vcf_to_merge=expand(VCF_DIR + "{sample}.dB.filt.recode.vcf.gz", sample=SAMPLES)
    output: VCF_MERGE_DIR + merged_sample_name + "_samples_to_merge.txt"
    run:
        with open(output[0], "w") as f:
            for sample in input.vcf_to_merge:
                f.write("{}\n".format(sample))

rule merge_vcf:
    """
    merge VCFs
    """
    input: VCF_MERGE_DIR + merged_sample_name + "_samples_to_merge.txt"
    output: VCF_MERGE_DIR + merged_sample_name + ".merge.vcf.gz"
    threads: 16
    shell:  """
        module purge
        module load bcftools
        bcftools merge --threads {threads} -l {input} -o {output} -O b
    """

rule extract_variants:
    """
    extract snps, and indels to a variant only VCF
    """
    input:  VCF_MERGE_DIR + merged_sample_name + ".merge.vcf.gz"
    output: VCF_MERGE_DIR + merged_sample_name + ".merge.variants.vcf"
    shell:  """
        module purge
        module load bcftools
        bcftools view -v snps,indels {input} -o {output}
    """

rule parse_vcf_to_multifasta:
    """
    parse snpEff annotated VCFs
    """
    input: VCF_MERGE_DIR + merged_sample_name + ".merge.variants.vcf"
    output: multi_fa=VCF_MERGE_DIR + merged_sample_name + ".merge.snps.indels.aln.fa"
    params: scripts=config["scripts"]
    threads: 16
    shell:  """
        module purge
        module load anaconda3
        python {params.scripts}/multiVCF_to_multifasta_v0.3.py -i {input} -o {output.multi_fa}
    """


# # Following version of rules are deprecated
# rule decompose_blocksub:
#     """
#     decompose blocksub, decomposes mnps into individual snps
#     """
#     input:  VCF_MERGE_DIR + merged_sample_name + ".merge.variants.vcf"
#     output: VCF_MERGE_DIR + merged_sample_name + ".merge.variants.Db.vcf"
#     shell: """
#     module purge
#     module load vt
#     vt decompose_blocksub -o {output} {input}
#     """

# rule normalize:
#     input:  VCF_MERGE_DIR + merged_sample_name + ".merge.variants.vcf"
#     output: VCF_MERGE_DIR + merged_sample_name + ".merge.variants.N.vcf"
#     params: reference=config["idx_bwa"]
#     shell: """
#     module purge
#     module load vt
#     vt normalize -r {params.reference} -o {output} {input}
#     """

# rule normalize:
#     input:  VCF_MERGE_DIR + merged_sample_name + ".merge.variants.vcf"
#     output: VCF_MERGE_DIR + merged_sample_name + ".merge.variants.N.vcf"
#     shell: """
#     module purge
#     module load bcftools
#     bcftools norm -m -any -c x {input} -o {output} 
#     """

# rule extract_snps:
#     input:  VCF_MERGE_DIR + merged_sample_name + ".merge.variants.N.vcf"
#     output: VCF_MERGE_DIR + merged_sample_name + ".merge.variants.N.snps.vcf"
#     shell:  """
#         module purge
#         module load bcftools
#         bcftools view -v snps {input} -o {output}
#     """
