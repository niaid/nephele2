
# calling specific helper functions
def get_vartype_arg(wildcards):
    return "--select-type-to-include {}".format(
        "SNP" if wildcards.vartype == "snvs" else "INDEL")

def get_filter_info(wildcards):  
    return [config["filtering"]["hard"][wildcards.vartype]["info"]]

def get_filter_genotype(wildcards):  
    return [config["filtering"]["hard"][wildcards.vartype]["genotype"]]

def get_filt_samples(wildcards):
    return expand("variants/filtered/{sample}.{vartype}.filt.vcf.gz",
                    vartype=["snvs", "indels"], **wildcards)

def get_filt_merge(wildcards):
    return " ".join("INPUT={}".format(f) for f in expand("variants/filtered/{sample}.{vartype}.filt.vcf.gz",
                    vartype=["snvs", "indels"], **wildcards))

def get_all_vcfs(wildcards):
    return expand("variants/filtered/{sample}.filt.GTset.vcf.gz", sample=samples.index)

# for just calling haplotype without using GVCF and combining
rule call_variants:
    input:
        bam=get_sample_bams,
        bai=get_sample_bai,
        ref=config["ref"]["genome"]
    output:
        vcf="variants/called_ind/{sample}.vcf.gz"
    params:
        ploidy=config["params"]["gatk"]["HaplotypeCaller"]["ploidy"]
    log:
        "log/gatk/haplotypecaller/{sample}.log"
    shell:"""
        module load gatk
        $EBROOTGATK/gatk --java-options "-Xmx4g" HaplotypeCaller  \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.vcf} \
            -ploidy {params.ploidy}
            """

# filtering on individual
rule select_calls:
    input:
        ref=config["ref"]["genome"],
        vcf="variants/called_ind/{sample}.vcf.gz"
    output:
        vcf=temp("variants/filtered/{sample}.{vartype}.vcf.gz"),
        tbi=temp("variants/filtered/{sample}.{vartype}.vcf.gz.tbi")
    params:
        extra=get_vartype_arg
    log:
        "log/gatk/selectvariants/{sample}.{vartype}.log"
    shell:"""
        module load gatk
        $EBROOTGATK/gatk SelectVariants -R {input.ref} -V {input.vcf} {params.extra} -O {output.vcf}
    """

# hard filter
rule hard_filter_calls:
    input:
        ref=config["ref"]["genome"],
        vcf="variants/filtered/{sample}.{vartype}.vcf.gz",
        tbi="variants/filtered/{sample}.{vartype}.vcf.gz.tbi"
    output:
        vcf=temp("variants/filtered/{sample}.{vartype}.filt.vcf.gz"),
        tbi=temp("variants/filtered/{sample}.{vartype}.filt.vcf.gz.tbi")
    params:
        filters_info=get_filter_info,
        filters_genotype=get_filter_genotype
    log:
        "log/gatk/variantfiltration/{sample}.{vartype}.log"
    shell:"""
        module load gatk
        $EBROOTGATK/gatk VariantFiltration \
           -R {input.ref} \
           -V {input.vcf} \
           -O {output.vcf} \
           --filter-name "snv-hard-filter" \
           --filter-expression "{params.filters_info}" \
           --genotype-filter-name "depth_filter" \
           --genotype-filter-expression "{params.filters_genotype}"
        """

# merge individual filtered
rule merge_ind_calls:
    input:
        vcf=get_filt_samples
    output:
        vcf="variants/filtered/{sample}.filt.vcf.gz"
    params:
        inputs=get_filt_merge
    log:
        "log/picard/merge-filtered.{sample}.log"
    shell:"""
        module load picard
        java -jar ${{EBROOTPICARD}}/picard.jar MergeVcfs \
            {params.inputs} \
            O={output}
    """

rule set_failed_GT:
    input:
        vcf="variants/filtered/{sample}.filt.vcf.gz"
    output:
        vcf="variants/filtered/{sample}.filt.GTset.vcf.gz"
    shell:"""
        module load bcftools
        bcftools filter -i 'FILTER=="PASS"' {input.vcf} -S . -o {output.vcf} -O b
        bcftools index {output.vcf} 
    """

# merge individual filtered
rule merge_vcfs:
    input:
        vcfs=expand("variants/filtered/{sample}.filt.GTset.vcf.gz", sample=samples.index)
    output:
        vcf="variants/combined/all.filt.GTset.vcf.gz"
    params:
        vcfs_each=get_all_vcfs
    log:
        "log/gatk/combinegvcfs.log"
    threads: 16
    shell:"""
        module load bcftools
        bcftools merge --threads {threads} -F x {params.vcfs_each} -o {output} -O b
    """

# parse PASSED from merge
rule merge_pass:
    input:
        vcf="variants/combined/all.filt.GTset.vcf.gz"
    output:
        vcf="variants/combined/all.filt.PASS.GTset.vcf.gz"
    log:
        "log/gatk/combinegvcfs.log"
    shell:"""
        module load bcftools
        bcftools view -f PASS {input.vcf} -o {output.vcf} -O b
    """

