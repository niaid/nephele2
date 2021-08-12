# rule vcf_to_tsv:
#     input:
#         "variants/filtered/all.p1.genotype.{resolution}.vars.hardfiltered.vcf.gz"
#     output:
#         report("tables/calls.{resolution}.tsv.gz", caption="../report/calls.{resolution}.rst", category="Calls")
#     conda:
#         "../envs/rbt.yaml"
#     shell:
#         "bcftools view --apply-filters PASS --output-type u {input} | "
#         "rbt vcf-to-txt -g --fmt DP AD --info ANN | "
#         "gzip > {output}"

rule vcf_to_tsv:
    input:
        "variants/filtered/all.p1.genotype.{resolution}.vars.hardfiltered.vcf.gz"
    output:
        report("variants/tables/calls.{resolution}.tsv.gz", caption="../report/calls.{resolution}.rst", category="Calls")
    conda:
        "../envs/rbt.yaml"
    shell:
        "bcftools view --output-type u {input} | "
        "rbt vcf-to-txt -g --fmt DP AD --info ANN | "
        "gzip > {output}"


# rule plot_stats:
#     input:
#         "tables/calls.tsv.gz"
#     output:
#         depths=report("plots/depths.svg", caption="../report/depths.rst", category="Plots"),
#         freqs=report("plots/allele-freqs.svg", caption="../report/freqs.rst", category="Plots")
#     conda:
#         "../envs/stats.yaml"
#     script:
#         "../scripts/plot-depths.py"