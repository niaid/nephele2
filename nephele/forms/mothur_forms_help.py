from nephele2 import tfvars

mothur_max_seq_length_desc = """
    Maximum length of sequences after reads are joined. If you leave the default value of 0, we will infer the correct maximum length from your data. If you wish to override the inferred value, you can set this field to a non-zero value.<br>
    Note: can cause job failure if set too low.
    """

mothur_difference_rank_desc = """
    The threshold of difference rank between rare and more abundant sequences.
    """

mothur_remove_lineage_desc = """
    Remove sequences that belong to taxa of certain lineages.
    """

mothur_remove_rare_otus_desc = """
    Remove OTUs at a specified rarity (number of observations).
    """

mothur_optimize_criteria_desc = """
    The <strong>optimize</strong> and <strong>criteria</strong> settings allow you set the <strong>start</strong>; and/or the <strong>end</strong> of your aligned sequences; as default, it would remove any sequence that starts after the position that 90% of the sequences, or ends before the position that 90% of the sequences. This ensures that the sequences align to the same region.
    """

mothur_reference_db_desc = """
    <strong>SILVA</strong>: the SILVA database will be used for taxonomic assignment.</br>
    <strong>HOMD</strong>: the HOMD database will be used for taxonomic assignment. Note: HOMD (Human Oral Microbiome Database) is a custom database for the analysis of human oral microbiome only.
    """

mothur_picrust_desc = """
    Checking this box to run <a href='http://picrust.github.io/picrust/' target="_blank" rel="noopener noreferrer">PICRUSt</a>.
    PICRUSt predicts metagenome functional content from marker genes. The analysis will
    produce two functional bar charts at taxa level 2 (Phylum) and 3 (Class).<br>
    Note: PICRUSt runs separately using Closed Reference OTU picking with the Greengenes 99 database.
    """

sampling_depth_desc = "The number of counts for filtering and subsampling the OTU table for downstream analysis. \
    Samples which have counts below this value will be removed from the downstream analysis. \
    The counts of the remaining samples will be subsampled to this value. If not specified, \
    it will be calculated automatically. \
    <a href=\""+tfvars.SERVER_ADDR+"/faq/#collapse12\" target=\"_blank\">See FAQ: Why are only some of my samples appearing in the downstream analysis for the amplicon pipelines?</a>"
