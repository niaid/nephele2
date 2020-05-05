"""
Defines the help text for the DADA2 pipeline form.
"""
from nephele2 import tfvars

ion_torrent_seq = """Checking this option sets the denoising parameters according to <a href="https://benjjneb.github.io/dada2/faq.html#can-i-use-dada2-with-my-454-or-ion-torrent-data" target="_blank" rel="noopener noreferrer">DADA2's suggested values for Ion Torrent data</a>. They also suggest the <code>trim left</code> parameter (below) be increased by 15 bp on top of any primer lengths.<br /><br />This option is in <span style="color:#FF0000">beta</span>, and has not been extensively tested. If you have Ion Torrent data, we are interested in your feedback - please <a href="mailto:nephelesupport@nih@@gov?Subject=Ion Torrent data" onmouseover="this.href=this.href.replace('@@','.')" onclick="this.href=this.href.replace('@@','.')" target="_top">email us</a>!"""

remove_chimera_desc = """
    Remove chimeric sequences. If primers are not trimmed
    (either prior to submission or using the left trim option), then we suggest unchecking this option.
    """

tax_method = """
    Method to be used for taxonomic assignment, either rdp or <a href=\"http://www2.decipher.codes/ClassificationFAQ.html\" target=\"_blank\">IDTAXA</a>.
    """

trim_left_fwd_desc = "The number of nucleotides to remove from the start of each read, forward and reverse. \
    The values should be chosen based on the lengths of primers used for sequencing. \
    If your data are untrimmed, this parameter is very important for the DADA2 pipeline. \
    <a href=\""+tfvars.SERVER_ADDR+"/faq/#collapse13\" target=\"_blank\">See this FAQ.</a>"

trim_left_rev_desc = "The number of nucleotides to remove from the start of each read, forward and reverse. \
    The values should be chosen based on the lengths of primers used for sequencing. \
    If your data are untrimmed, this parameter is very important for the DADA2 pipeline. \
    <a href=\""+tfvars.SERVER_ADDR+"/faq/#collapse13\" target=\"_blank\">See this FAQ.</a>"

trim_overhang_desc = """
    After merging paired end reads, trim sequence which overhangs the start of each read.
    If amplicons are shorter than read length, e.g. 16S V4 region, we suggest checking this option.
    """

just_concatenate= """
    Concatenate paired reads instead of merging.
    """

tuncq_desc = """
    Truncate reads at the first instance of a quality score less than or equal to this value.
    """

truncLen_fwd_desc = """
    The length at which to truncate reads, forward and reverse. Reads shorter than these lengths are discarded. If set to 0, reads are not truncated.
    """

truncLen_rev_desc = """
    The length at which to truncate reads, forward and reverse. Reads shorter than these lengths are discarded. If set to 0, reads are not truncated.
    """

max_mismatches_desc = """
    The maximum number of mismatches allowed in the overlap region when merging read pairs.
"""

maxee_desc = """
    After truncation, reads with higher than this many “expected errors” will be discarded. Expected errors are calculated from the nominal definition of the quality score: EE = sum(10^(-Q/10)).
"""

ref_db_desc = """
    <strong>SILVA</strong>: the SILVA database will be used for taxonomic assignment.</br>
    <strong>HOMD</strong>: the HOMD database will be used for taxonomic assignment. Note: HOMD (Human Oral Microbiome Database) is a custom database for the analysis of human oral microbiome only.
    IDTAXA will use its own SILVA v132 database.
"""

sampling_depth_desc = "The number of counts for filtering and subsampling the OTU table for downstream analysis. \
    Samples which have counts below this value will be removed from the downstream analysis. \
    The counts of the remaining samples will be subsampled to this value. If not specified, \
    it will be calculated automatically. \
    <a href=\""+tfvars.SERVER_ADDR+"/faq/#collapse12\" target=\"_blank\">See FAQ: Why are only some of my samples appearing in the downstream analysis for the amplicon pipelines?</a>"
