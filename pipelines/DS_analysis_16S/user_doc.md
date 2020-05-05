---
panflute-filters: [addlinkattr, inserttoc]
panflute-path: '../../resources/Rdocs'
codebox: false #do we need syntax highlighting for code?
output:
  html_document:
    mathjax: null
    section_divs: no
    smart: no
    md_extensions: "-smart-autolink_bare_uris"
    highlight: null
    theme: null
    template: ../../resources/Rdocs/user_doc_template.html #change path to wherever it is
    toc_depth: 3
    toc: yes
    self_contained: no
---

# Amplicon Downstream Analysis Pipeline {#title}

The Downstream Analysis pipeline provides sample observation and taxonomic summaries and diversity analyses of an OTU table using [QIIME 2](https://qiime2.org).

## Input Files and Parameters

- **BIOM File:** The biom file contains the OTU and taxonomy tables to be analyzed.  This pipeline accepts [BIOM V1](http://biom-format.org/documentation/format_versions/biom-1.0.html) or QIIME's BIOM V2 format.  All biom files produced by Nephele's amplicon pipelines should work.
- **Mapping File:** The mapping file contains the metadata about the samples which will be used in the analysis.  The mapping file format is the same as that used by the amplicon pipelines (<a href="{{ se_map_url }}" target="_blank" rel="noopener noreferrer">SE</a> or <a href="{{ pe_map_url }}" target="_blank" rel="noopener noreferrer">PE templates</a>).  The FASTQ file columns will be ignored.  Only samples listed in the mapping file will be included in the pipeline results.  So, if you want to analyze only a subset of the samples in your biom file, you can submit the mapping file with only those samples, and the others will be excluded from the analysis.
- **Sampling depth:**  The total frequency that each sample should be rarefied to prior to computing diversity metrics.  Samples with counts below this value will be removed from the analysis.  If you are using the output of a Nephele amplicon pipeline, you may find it useful to consult *otu_summary_table.txt* which lists the sample counts.
- **Alpha group significance:** Run alpha diversity statistical comparisons between groups, and produce alpha diversity plots.  This step runs after samples with counts below your sampling depth are filtered from your biom and mapping files.  It requires the filtered mapping file's Treatment Group column to contain at least 2 groups, and each group to contain at least 2 samples.

  If you are not sure which samples and groups will remain after the sampling depth filtering, run the job first with this option unchecked, and review the *summary.qzv* file.  If you are using the output of a Nephele amplicon pipeline, you may find it useful to consult *otu_summary_table.txt* which lists the sample counts.

## QIIME 2

Nephele uses the [QIIME 2 v2018.6](https://docs.qiime2.org/2018.6/) Artifact [API](https://docs.qiime2.org/2018.6/interfaces/artifact-api/).

## Pipeline Steps & Output files {#output-files}

The output `.qza` and `.qzv` files can be viewed on QIIME 2's [view page](https://view.qiime2.org). See QIIME 2's information about [the output formats](https://docs.qiime2.org/2018.6/concepts/#data-files-visualizations) and for [help](https://view.qiime2.org/about) with the view page.  Where possible, the `.qza` artifacts are also exported to their native format in directories of the same name.  The following plugin methods are used:

1. [**feature-table summarize**](https://docs.qiime2.org/2018.6/plugins/available/feature-table/summarize/) which produces a summary of the counts along with the sample metadata.
   - *summary.qzv*: summary visualization
2. [**taxa barplot**](https://docs.qiime2.org/2018.6/plugins/available/taxa/barplot/) which produces barplots of the taxonomies of the samples with counts above the sampling depth.  [filter-samples](https://docs.qiime2.org/2018.6/plugins/available/feature-table/filter-samples/) is used to filter these samples.
   - *barplot.qzv*
3. [**diversity core-metrics**](https://docs.qiime2.org/2018.6/plugins/available/diversity/core-metrics/) which applies a set of non-phylogenetic diversity metrics to the count data after rarefaction.
   - *rarefied_table.qza, rarefied_table/feature-table.biom*: rarefied feature table
   - *observed_otus_vector.qza*: vector of observed OTU/ASV values by sample
   - *shannon_vector.qza*: vector of Shannon diversity values by sample
   - *evenness_vector.qza*: vector of Pielou's evenness values by sample
   - *jaccard_distance_matrix.qza,  jaccard_pcoa_results.qza,  jaccard_emperor.qzv*: Jaccard distance matrix and coordinates for the resulting PCoA ordination.  The graph is plotted using [Emperor](https://biocore.github.io/emperor/).
   - *bray_curtis_distance_matrix.qza, bray_curtis_pcoa_results.qza, bray_curtis_emperor.qzv*: Bray-Curtis distance matrix and coordinates for the resulting PCoA ordination.  The graph is plotted using [Emperor](https://biocore.github.io/emperor/).
4. [**diversity alpha-group-significance**](https://docs.qiime2.org/2018.6/plugins/available/diversity/alpha-group-significance/) which does statistical comparisons of the alpha diversity indexes based on the sample metadata groups.
   - *alpha_group_significance_evenness.qzv*:  comparison of Pielou (evenness) index between groups
   - *alpha_group_significance_shannon.qzv*:  comparison of Shannon index between groups

## Tools and References
- Bolyen, Evan, et al. "Reproducible, Interactive, Scalable and Extensible Microbiome Data Science Using QIIME 2." *Nature Biotechnology*, July 2019, doi:[10.1038/s41587-019-0209-9](https://doi.org/10.1038/s41587-019-0209-9).
- McDonald, Daniel, et al. "The Biological Observation Matrix (BIOM) Format or: How I Learned to Stop Worrying and Love the Ome-Ome." *GigaScience*, vol. 1, July 2012, p. 7. *BioMed Central*, doi:[10.1186/2047-217X-1-7](https://doi.org/10.1186/2047-217X-1-7).
- Kruskal, William H., and W. Allen Wallis. "Use of Ranks in One-Criterion Variance Analysis." *Journal of the American Statistical Association*, vol. 47, no. 260, 1952, pp. 583-621. *JSTOR*, doi:[10.2307/2280779](https://doi.org/10.2307/2280779).
- Vázquez-Baeza, Yoshiki, Antonio Gonzalez, et al. "Bringing the Dynamic Microbiome to Life with Animations." *Cell Host &amp; Microbe*, vol. 21, no. 1, Jan. 2017, pp. 7-10. *ScienceDirect*, doi:[10.1016/j.chom.2016.12.009](https://doi.org/10.1016/j.chom.2016.12.009).
- Vázquez-Baeza, Yoshiki, Meg Pirrung, et al. "EMPeror: A Tool for Visualizing High-Throughput Microbial Community Data." *GigaScience*, vol. 2, no. 1, Dec. 2013. *academic.oup.com*, doi:[10.1186/2047-217X-2-16](https://doi.org/10.1186/2047-217X-2-16).
- Weiss, Sophie, et al. "Normalization and Microbial Differential Abundance Strategies Depend upon Data Characteristics." *Microbiome*, vol. 5, no. 1, Mar. 2017, p. 27. *Crossref*, doi:[10.1186/s40168-017-0237-y](https://doi.org/10.1186/s40168-017-0237-y).
