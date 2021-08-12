User docs
================

-   [ampvis2 and Plotly](#ampvis2-and-plotly)
-   [Plots](#plots)
-   [Output Files](#output-files)
-   [Tools and References](#tools-and-references)

The 16S visualization pipeline runs automatically at the end of Nephele’s 16S amplicon pipelines and produces some basic plots for rarefaction, alpha and beta diversity and heatmaps of taxonomic composition. For additional visualizations and more fine-grained diversity analysis, please see our separate <a class="" href="{{ url_for('show_pipes_guide', _anchor='da_pipes') }}">downstream analysis (*DA*) pipeline</a>.

<div class="mt-5 mb-5">

## ampvis2 and Plotly

Nephele uses the [ampvis2 R package](https://madsalbertsen.github.io/ampvis2/) for statistical computation. We also make use of the [Plotly R interface](https://plot.ly/r/) for the [plotly.js](https://plot.ly) charting library. The plots can be minimally edited interactively. However, they also have a link to export the data to the [Plotly Chart Studio](https://plotly.com/chart-studio/) which allows for [making a variety of charts](https://plotly.com/chart-studio-help/tutorials/). We have a tutorial for <a href="{{ url_for('show_tutorials') }}">making simple edits to the plots here</a>.

</div>

<div class="mt-5 mb-5">

## Plots

### Rarefaction Curve

The rarefaction curves are made with the [amp_rarecurve](https://madsalbertsen.github.io/ampvis2/reference/amp_rarecurve.html) function. The table is written to *rarecurve.txt* and the plot to *rarecurve.html*.

``` r
rarecurve <- amp_rarecurve(amp, color_by = "TreatmentGroup")
```

### Filtering and rarefying/subsampling OTU table

The **`sampling depth`** is used to filter out samples with few counts. Samples with counts which fall below the `sampling depth` are removed from the OTU table using [amp_subset_samples](https://madsalbertsen.github.io/ampvis2/reference/amp_subset_samples.html). The names of samples that are removed are output to *samples_being_ignored.txt*.

Additionally, the `sampling depth` is used to rarefy the counts for the Bray-Curtis PCoA and the alpha diversity plots. The OTU table is rarefied using [rrarefy](https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/rarefy) from the vegan R package. The table of the filtered, rarefied counts is saved as *rarefied_OTU_table.txt*.

If you do not provide a `sampling depth`, the default value is **10000**.

``` r
ampsub <- amp_subset_samples(amp, minreads = sampdepth)
otu <- rrarefy(t(ampsub$abund), sampdepth)
amprare$abund <- t(otu)
```

### Heatmap

The interactive heatmap is implemented using the [morpheus R API](https://github.com/cmap/morpheus.R) developed at the Broad Institute. [Documentation](https://software.broadinstitute.org/morpheus/documentation.html) for how to use the heatmap can be found on the [morpheus website](https://software.broadinstitute.org/morpheus/).

The heatmap, *seq_heatmap.html*, is made from the raw OTU table with the counts normalized to 100 to represent the relative abundances using [amp_subset_samples](https://madsalbertsen.github.io/ampvis2/reference/amp_subset_samples.html) before the heatmap is made with morpheus. If there are too many OTUs or sequence variants, then the heatmap is made at the species level instead. A heatmap is also made from the rarefied counts, *seq_heatmap_rarefied.html*.

``` r
amptax <- amp_subset_samples(amp, normalise = TRUE)
heatmap <- morpheus(amp$abund, columns = columns, rows = rows,
    columnColorModel = list(type = as.list(colors)),
    colorScheme = list(scalingMode = "fixed", stepped = FALSE),
    columnAnnotations = amptax$metadata, rowAnnotations = amptax$tax)
```

### PCoA

Principal coordinates analysis using binomial and Bray-Curtis distances is carried out using [amp_ordinate](https://madsalbertsen.github.io/ampvis2/reference/amp_ordinate.html). For more information on the formulae for the distance measures, see [`vegdist`](https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/vegdist).

The binomial distance is able to handle varying sample sizes, so the raw counts from the OTU table are used. For the Bray-Curtis distance, the rarefied counts are used. The coordinates from the plots are written to *pcoa_binomial.txt* and *pcoa_bray.txt* and the plots to *pcoa_binomial.html* and *pcoa_bray.html*. At least 3 samples are needed for these plots.

``` r
pcoa_binomial <- amp_ordinate(amp, filter_species = 0.01, type = "PCOA",
    distmeasure = "binomial", sample_color_by = "TreatmentGroup",
    detailed_output = TRUE, transform = "none")

pcoa_bray <- amp_ordinate(amprare, filter_species = 0.01, type = "PCOA",
    distmeasure = "bray", sample_color_by = "TreatmentGroup",
    detailed_output = TRUE, transform = "none")
```

### Alpha diversity

The Shannon diversity and Chao species richness are computed using [amp_alphadiv](https://madsalbertsen.github.io/ampvis2/reference/amp_alphadiv.html). The rarefied counts are used for this computation. The diversity values are saved in *alphadiv.txt*, and boxplots are output to *alphadiv.html*. At least 3 samples are needed to produce these plots.

``` r
alphadiv <- amp_alphadiv(amprare, measure = "shannon", richness = TRUE, rarefy = sampdepth)
```

</div>

<div class="mt-5 mb-5">

## Output Files

Complete descriptions of the output files can be found in the [Plots section above](#plots). To learn how to edit the plots, see the <a href="{{ url_for('show_tutorials') }}">visualization tutorial</a>.

-   *rarecurve.html*: rarefaction curve plot
-   *rarecurve.txt*: tabular data used to make the rarefaction curve plot
-   *seq_heatmap\*.html*: heatmap of OTU/sequence variant abundances
-   *samples_being_ignored.txt*: list of samples removed from the analysis
-   *pcoa\_\*.html*: PCoA plots
-   *pcoa\_\*.txt*: tabular data used to make the PCoA plots
-   *rarefied_OTU_table.txt*: rarefied OTU table used for Bray-Curtis PCoA and alpha diversity plots
-   *alphadiv.html*: alpha diversity boxplots
-   *alphadiv.txt*: tabular data used to make the alpha diversity plots

</div>

<div class="mt-5 mb-5">

## Tools and References

Andersen KS, Kirkegaard RH, Karst SM, Albertsen M (2018). “ampvis2: an R package toand visualise 16S rRNA amplicon data.” *bioRxiv*. <https://doi.org/10.1101/299537>.
<p>
Gould J (2019). <em>morpheus: Interactive heat maps using ‘morpheus.js’ and ‘htmlwidgets’</em>. R package version 0.1.1.1, <a href="https://github.com/cmap/morpheus.R">https://github.com/cmap/morpheus.R</a>.
</p>
<p>
McMurdie PJ and Paulson JN (2016). <em>biomformat: An interface package for the BIOM file format</em>. <a href="https://github.com/joey711/biomformat/" target="_blank" rel="noopener noreferrer">https://github.com/joey711/biomformat/</a>.
</p>
<p>
Oksanen J, Blanchet FG, Friendly M, Kindt R, Legendre P, McGlinn D, Minchin PR, O’Hara RB, Simpson GL, Solymos P, Stevens MHH, Szoecs E, Wagner H (2019). <em>vegan: Community Ecology Package</em>. R package version 2.5-4, <a href="https://CRAN.R-project.org/package=vegan" target="_blank" rel="noopener noreferrer">https://CRAN.R-project.org/package=vegan</a>.
</p>

Sievert C (2020). *Interactive Web-Based Data Visualization with R, plotly, and shiny*. Chapman and Hall/CRC. ISBN 9781138331457, <https://plotly-r.com>.

</div>
