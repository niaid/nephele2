
rule trim_reads_pe:
    input:
        unpack(get_fastq_gz)
    output:
        r1="data/trimmed/{sample}.trim.1.fastq.gz",
        r2="data/trimmed/{sample}.trim.2.fastq.gz",
        r1_unpaired="data/trimmed/{sample}.1.unpaired.fastq.gz",
        r2_unpaired="data/trimmed/{sample}.2.unpaired.fastq.gz",
        trimlog="data/trimmed/{sample}.trimlog.txt"
    params:
        extra=lambda w, output: "-trimlog {}".format(output.trimlog),
        trimmer=" ".join(config["params"]["trimmomatic"]["pe"]["trimmer"]),
        adapter=config["params"]["trimmomatic"]["illumina_clip"]
    threads: 16
    shell:"""
        module purge
        module load trimmomatic
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads {threads} {params.extra} {input.r1} {input.r2} {output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired} ILLUMINACLIP:{params.adapter} {params.trimmer} 
    """

rule trim_reads_se:
    input:
        unpack(get_fastq_gz)
    output:
        r1="data/trimmed/{sample}.trim.fastq.gz",
        trimlog="data/trimmed/{sample}.trimlog.txt"
    params:
        extra=lambda w, output: "-trimlog {}".format(output.trimlog),
        trimmer=" ".join(config["params"]["trimmomatic"]["se"]["trimmer"]),
        adapter=config["params"]["trimmomatic"]["illumina_clip"]
    threads: 16
    shell:"""
        module purge
        module load trimmomatic
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar SE -threads {threads} {params.extra} {input.r1} {output.r1} ILLUMINACLIP:{params.adapter} {params.trimmer} 
    """

rule map_reads:
    input:
        reads=get_trimmed_reads
    output:
        temp("mapped/{sample}.sort.bam")
    log:
        "log/bwa_mem/{sample}.log"
    params:
        index=config["ref"]["genome"],
        extra=get_read_group
    threads: 8

    shell:"""
        echo {input.reads}
        module load bwa samtools
        bwa mem -t {threads} {params.extra} {params.index} {input.reads} | samtools sort - > {output} 2> {log}
        """

rule mark_duplicates:
    input:  
        "mapped/{sample}.sort.bam"   
    output: 
        bam="mapped/{sample}.sort.dedup.bam"
    shell:"""
        module load samtools
        samtools rmdup {input} {output.bam}
        samtools index {output.bam}
        """

rule per_base_coverage:
    input:
        "mapped/{sample}.sort.dedup.bam"
    output: 
        "qc/coverage/{sample}.full.coverage.txt"
    params:
        index=config["ref"]["bedtools_genome"]
    shell:"""
        module load bedtools
        bedtools genomecov -ibam {input} -g {params.index} -d > {output}
        """

rule coverage_QC:
    input:
        coverage_files=expand("qc/coverage/{sample}.full.coverage.txt", sample=samples.index),
        sample_table=config["samples"]
    output:
        coverage_QC_file=config["samples_QC"],
        where_at="log/coverage_sample_check.txt"
    script:
        "../scripts/coverage_process.py"
        
