shell.prefix("set -eo pipefail; ")

##### Wildcard constraints #####
wildcard_constraints:
    vartype="snvs|indels",
    sample="|".join(samples.index)

##### Helper functions #####

def get_fastq(wildcards):
    """
    Get fastq files of given sample
    """
    fastqs = samples.loc[wildcards.sample, ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}

def get_fastq_gz(wildcards):
    """
    Get fastq.gz  files of given sample.
    """
    fastqs = samples.loc[wildcards.sample, ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1 + ".gz", "r2": fastqs.fq2 + ".gz"}
    return {"r1": fastqs.fq1 + ".gz"}

def get_read_group(wildcards):
    """
    sample name and platform in read group
    """
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform=samples.loc[wildcards.sample, "platform"])

def is_single_end(sample):
    """
    Return True if sample is single end
    """
    return pd.isnull(samples.loc[sample, "fq2"])

def get_trimmed_reads(wildcards):
    """
    Get trimmed reads of given sample
    """
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("data/trimmed/{sample}.trim.{group}.fastq.gz", 
            group=[1, 2], **wildcards)

    # single end sample
    return "data/trimmed/{sample}.trim.fastq.gz".format(**wildcards)

def get_sample_bams(wildcards):
    """
    Get all aligned reads of given sample.
    """
    return expand("mapped/{sample}.sort.dedup.bam",
                  sample=wildcards.sample)

def get_sample_bai(wildcards):
    """
    Get all indexed bam.
    """
    return expand("mapped/{sample}.sort.dedup.bam.bai",
                  sample=wildcards.sample)

def get_merge_file_prefix():
    return config["samples_QC"].rsplit(".", 1)[0]
    
