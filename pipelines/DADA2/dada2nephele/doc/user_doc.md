DADA2 Pipeline
================

  - [R packages](#r-packages)
  - [User Options](#user-options)
  - [Databases](#databases)
  - [Pipeline steps](#pipeline-steps)
  - [Output Files](#output-files)
  - [Tools and References](#tools-and-references)

## R packages

Nephele runs the [DADA2 R
package](https://benjjneb.github.io/dada2/index.html) v1.10 following
the steps in the package authors' [Big Data
workflow](https://benjjneb.github.io/dada2/bigdata.html) including
optional use of [DECIPHER](http://www2.decipher.codes) package v2.10. We
make some minor modifications of the parameters used. Our pipeline is
[outlined below](#pipeline-steps). **If you are new to DADA2, it might
be helpful to read through the [DADA2
Tutorial](https://benjjneb.github.io/dada2/tutorial.html).**

## User Options

  - **Ion Torrent Data -** <strong style="color:#FF0000">Beta</strong>:
    By default, DADA2 is trained to work on Illumina data. Checking this
    option sets the denoising parameters according to [DADA2's suggested
    values for Ion Torrent
    data](https://benjjneb.github.io/dada2/faq.html#can-i-use-dada2-with-my-454-or-ion-torrent-data).
    They also suggest the `trim left` parameter be increased by 15 bp
    (on top of any primer lengths). This option is in
    <span style="color:#FF0000">beta</span>, and has not been
    extensively tested. If you have Ion Torrent data, we are interested
    in your feedback - please
    <a href="mailto:nephelesupport@nih@@gov?Subject=Ion Torrent data" onmouseover="this.href=this.href.replace('@@','.')" onclick="this.href=this.href.replace('@@','.')" target="_top">email
    us</a>\!

#### Filter and Trim

  - **Trim left**: The number of nucleotides to remove from the start of
    each read, forward and reverse. The values should be chosen based on
    the lengths of primers used for sequencing. If your data are
    untrimmed, this parameter is very important for the DADA2 pipeline.
    See this
    <a href="{{ url_for('show_FAQ', _anchor='collapse13') }}">FAQ</a>
    (Default: 20).
  - **Truncation quality score**: Truncate reads at the first instance
    of a quality score less than or equal to this value. (Default: 4).
  - **Truncation length**: The length at which to truncate reads,
    forward and reverse. Reads shorter than these lengths are discarded.
    If set to 0, reads are not truncated. If both trim left and
    truncation length are set, the filtered reads will have length =
    *truncation length - trim left*. (Default: 0).
  - **Maximum expected errors (maxEE)**: After truncation, reads with
    higher than this many "expected errors" will be discarded. Expected
    errors are calculated from the nominal definition of the quality
    score: EE = sum(10^(-Q/10)). (Default: 5).

#### Merge Pairs

For *paired-end data only*.

  - **Just concatenate**: Concatenate paired reads instead of merging.
    (Default: FALSE)
  - **Maximum mismatches**: The maximum number of mismatches allowed in
    the overlap region when merging read pairs. (Default: 0).
  - **Trim overhanging sequence**: After merging paired end reads, trim
    sequence which overhangs the start of each read. If amplicons are
    shorter than read length, e.g. 16S V4 region, we suggest checking
    this option. (Logical. Default: FALSE).

#### Analysis

  - **Chimera removal**: Remove chimeric sequences. If primers are not
    trimmed (either prior to submission or using the trim left option),
    then we suggest unchecking this option. (Default: True).
  - **Taxonomic assigment**: Method to be used for taxonomic assignment,
    either `rdp` or
    [`IDTAXA`](http://www2.decipher.codes/ClassificationFAQ.html).
    (Default: rdp)
  - **Reference database**: Reference database to be used for taxonomic
    assignment. `IDTAXA` will use its own SILVA v132 database. See
    [Databases below](#databases).
  - **Sampling depth**: The number of counts for filtering and
    subsampling the OTU table for downstream analysis. Samples which
    have counts below this value will be removed from the downstream
    analysis. The counts of the remaining samples will be subsampled to
    this value. If not specified, it will be calculated automatically.
    (See <a href="{{ url_for('show_FAQ', _anchor='collapse18') }}">this
    FAQ</a>)

## Databases

  - [SILVA v132
    database](https://benjjneb.github.io/dada2/training.html)
  - [Human Oral Microbiome Database
    (HOMD)](http://www.homd.org/index.php) v15.1 formatted for DADA2
  - [Greengenes](https://benjjneb.github.io/dada2/training.html) v13.8
  - For `IDTAXA`, we use the authors' [modified SILVA v132 SSU trained
    classifier](http://www2.decipher.codes/Downloads.html). More
    information in [the DECIPHER
    FAQ](http://www2.decipher.codes/ClassificationFAQ.html).

## Pipeline steps

1)  [Plot quality
    profiles](https://rdrr.io/bioc/dada2/man/plotQualityProfile.html) of
    forward and reverse reads. These graphs are saved as
    *qualityProfile\_R1.pdf* and *qualityProfile\_R2.pdf*.
    
    ``` r
    pqp1 <- plotQualityProfile(file.path(datadir, r1))
    pqp2 <- plotQualityProfile(file.path(datadir, r2))
    ```

2)  Preprocess sequence data with
    [filterAndTrim](https://rdrr.io/bioc/dada2/man/filterAndTrim.html).
    The **maxEE**, **trimLeft**, **truncQ**, and **truncLen** parameters
    can be set by the user (defaults used below as example). The
    filtered sequence files, *\*\_trim.fastq.gz*, are output to the
    *filtered\_data* directory.
    
    ``` r
    filterAndTrim(fwd = file.path(datadir, readslist$R1), filt = file.path(filt.dir, trimlist$R1), 
        rev = file.path(datadir, readslist$R2), filt.rev = file.path(filt.dir, trimlist$R2), 
        maxEE = 5, trimLeft = 20, truncQ = 4, truncLen = 0, rm.phix = TRUE, compress = TRUE, 
        verbose = TRUE, multithread = nthread, minLen = 50)
    ```

3)  Learn the error rates with
    [learnErrors](https://rdrr.io/bioc/dada2/man/learnErrors.html). The
    **nbases** parameter is set to **1e+08**. The error rate graphs made
    with [plotErrors](https://rdrr.io/bioc/dada2/man/plotErrors.html)
    are saved as *errorRate\_R1.pdf*, *errorRate\_R2.pdf*. The error
    profiles, `err`, are also saved as a list R binary object in the
    *intermediate\_files* directory.
    
    ``` r
    errR1 <- learnErrors(r1, multithread = nthread, nbases = nbases, randomize = TRUE)
    pe1 <- plotErrors(errR1, nominalQ = TRUE)
    ```

4)  Dereplicate reads with
    [derepFastq](https://rdrr.io/bioc/dada2/man/derepFastq.html) and run
    the [dada](https://rdrr.io/bioc/dada2/man/dada.html)
    sequence-variant inference algorithm.
    
    ``` r
    derepR1 <- derepFastq(r1[sample], verbose = TRUE)
    ddR1 <- dada(derepR1, err = errR1, multithread = nthread, verbose = 1)
    ```

5)  For paired-end data, merge the overlapping denoised reads with
    [mergePairs](https://rdrr.io/bioc/dada2/man/mergePairs.html). The
    **minOverlap** parameter is set to **12**. **trimOverhang**,
    **justConcatenate**, and **maxMismatch** are set by the user. The
    sequence table, `seqtab`, containing the final amplicon sequence
    variants (ASVs), is saved as an R binary object to the
    *intermediate\_files* directory.
    
    ``` r
    mergePairs(dd$R1, derep$R1, dd$R2, derep$R2, verbose = TRUE, minOverlap = 12, trimOverhang = FALSE, 
        justConcatenate = FALSE, maxMismatch = 0)
    seqtab <- makeSequenceTable(mergedReads)
    ```

6)  Filter out ASVs of length less than 75 bp. The sequence table is
    saved as *seqtab\_min75.rds*. Also, filter out chimeras with
    [removeBimeraDenovo](https://rdrr.io/bioc/dada2/man/removeBimeraDenovo.html),
    if the option is chosen.
    
    ``` r
    seqtab <- seqtab[, which(seqlengths >= 75)]
    seqtabnochimera <- removeBimeraDenovo(seqtab, verbose = TRUE, multithread = nthread)
    ```

7)  Classify the remaining ASVs taxonomically with
    
      - **rdp using
        [assignTaxonomy](https://rdrr.io/bioc/dada2/man/assignTaxonomy.html)**
        (default). The **minBoot** parameter for minimum bootstrap
        confidence is set to **80** and **tryRC** is set to **TRUE**.
        This genus level result is saved as *taxa.biom*. Add species
        annotation to the taxonomic identification using
        [addSpecies](https://rdrr.io/bioc/dada2/man/addSpecies.html).
        This final result is saved as a biom file *taxa\_species.biom*.
        
        ``` r
        taxa <- assignTaxonomy(seqtab, refdb, multithread = nthread, minBoot = 80, tryRC = TRUE, 
            verbose = TRUE)
        taxa.species <- addSpecies(taxa, refdb_species, verbose = TRUE, tryRC = TRUE)
        ```
        
        For PE data, if the
        [mergePairs](https://rdrr.io/bioc/dada2/man/mergePairs.html)
        **justConcatenate** option is checked, species annotation will
        only be done [using the forward reads
        (R1)](https://github.com/benjjneb/dada2/issues/529#issuecomment-408171883).
    
      - **IDTAXA using
        [IdTaxa](https://rdrr.io/bioc/DECIPHER/man/IdTaxa.html)** from
        the DECIPHER R package. The final result will be saved as
        *taxa.biom*
        
        ``` r
        dna <- DNAStringSet(getSequences(seqtab))
        ids <- IdTaxa(dna, trainingSet, strand = "both", processors = nthread, verbose = T)
        ```

8)  The final results are also saved as a tab-separated text file
    *OTU\_table.txt*. The final sequence variants used for taxonomic
    classification are output as *seq.fasta*. A summary of the counts in
    the OTU table is saved to *otu\_summary\_table.txt*.

## Output Files

Complete descriptions of the intermediate and final output files can be
found in the [Pipeline Steps above](#pipeline-steps).

  - *OTU\_table.txt*: tab-separated text file of ASV counts and
    taxonomic assignment
  - *seq.fasta*: FASTA file of amplicon sequence variants
  - *taxa\_species.biom*: taxonomic assignment at the species level
    [BIOM V1
    format](http://biom-format.org/documentation/format_versions/biom-1.0.html)
  - *otu\_summary\_table.txt*: summary of the sequence variant counts by
    sample
  - *taxa.biom*: [BIOM
    V1](http://biom-format.org/documentation/format_versions/biom-1.0.html)
    file of sequence variant counts with taxonomic assignment at the
    genus level
  - *taxonomy\_table.txt*: tab-separated taxonomy file suitable for
    <a href="{{ url_for('show_qiime2_tutorial') }}">importing into
    QIIME2</a>
  - *errorRate\_R1/2.pdf*: error profile plots
  - *qualityProfile\_R1/2.pdf*: quality profile plots
  - **filtered\_data**: trimmed sequence files
  - **intermediate\_files**: intermediate files produced by the
    pipeline; useful for debugging
  - **graphs**: output of the
    <a href="{{ url_for('show_vis_details') }}">visualization
    pipeline</a>

## Tools and References

1.  Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA and Holmes
    SP (2016). "DADA2: High-resolution sample inference from Illumina
    amplicon data." <em>Nature Methods</em>, <b>13</b>, pp. 581-583.
    doi:
    <a href="http://doi.org/10.1038/nmeth.3869" target="_blank" rel="noopener noreferrer">10.1038/nmeth.3869</a>.

2.  McMurdie PJ and Paulson JN (2016). <em>biomformat: An interface
    package for the BIOM file format</em>.
    <a href="https://github.com/joey711/biomformat/" target="_blank" rel="noopener noreferrer">https://github.com/joey711/biomformat/</a>.

3.  Microsoft and Weston S (2017). <em>foreach: Provides Foreach Looping
    Construct for R</em>. R package version 1.4.4,
    <a href="https://CRAN.R-project.org/package=foreach" target="_blank" rel="noopener noreferrer">https://CRAN.R-project.org/package=foreach</a>.

4.  Murali, A., Bhargava, A., and Wright, E. S. (2018). IDTAXA: a novel
    approach for accurate taxonomic classification of microbiome
    sequences. <em>Microbiome</em>, <b>6</b>(1). doi:
    [10.1186/s40168-018-0521-5](https://doi.org/10.1186/s40168-018-0521-5).

5.  Quast C., Pruesse E., Yilmaz P., Gerken, J., Schweer T., Yarza P.,
    Peplies, J., Glöckner, F. O. (2013). "The SILVA ribosomal RNA gene
    database project: Improved data processing and web-based tools."
    *Nucleic Acids Research*, **41**(D1), D590-D596. doi:
    [10.1093/nar/gks1219](https://doi.org/10.1093/nar/gks1219).

6.  Escapa, I. F., Chen, T., Huang, Y., Gajare, P., Dewhirst, F. E., &
    Lemon, K. P. (2018). "New Insights into Human Nostril Microbiome
    from the Expanded Human Oral Microbiome Database (eHOMD): a Resource
    for the Microbiome of the Human Aerodigestive Tract."
    <em>MSystems</em>, <b>3</b>(6), e00187-18. doi:
    [10.1128/mSystems.00187-18](https://doi.org/10.1128/mSystems.00187-18).

7.  DeSantis, T. Z., et al. "Greengenes, a Chimera-Checked 16S rRNA Gene
    Database and Workbench Compatible with ARB." *Applied and
    Environmental Microbiology*, vol. 72, no. 7, July 2006, pp. 5069–72.
    *aem.asm.org*, doi:
    [10.1128/AEM.03006-05](https://doi.org/10.1128/AEM.03006-05).
