# de novo assembly

def paired_status(sample):
    """
    Return True if sample is single end
    """
    return samples.loc[sample, "library"] == "SINGLE"

rule pilon:
    """
    pilon assembly 
    """
    input:
        bam="mapped/{sample}.sort.dedup.bam",
        ref=config["ref"]["genome"]
    output: "pilon/{sample}.fasta"
    params: pilon="{sample}",
            pilon_dir="pilon",
            ref_name=config['ref']['name']
    threads: 16
    shell:"""
        module purge
        module load pilon/1.23
        java -Xmx16G -jar $EBROOTPILON/pilon-1.23.jar --genome {input.ref} --bam {input.bam} --output {params.pilon} --outdir {params.pilon_dir}
        sed -E -i "s/{params.ref_name}_pilon/{params.pilon}/g" {output}
    """

