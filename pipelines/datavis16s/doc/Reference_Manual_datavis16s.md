Package 'datavis16s'
================
August 03, 2021


```
Package: datavis16s
Title: Graphs for Nephele 16S Pipelines
Version: 0.1.3
Date: 2021-08-03
Authors@R (parsed):
    * Poorani Subramanian <poorani.subramanian@nih.gov> [aut, cre]
Description: betterbetterplots!
License: file LICENSE
URL:
    https://github.niaid.nih.gov/bcbb/nephele2/tree/master/pipelines/datavis16s
Depends:
    R (>= 3.4.0)
Imports:
    ampvis2,
    biomformat,
    data.table,
    ggplot2,
    htmltools,
    htmlwidgets,
    jsonlite,
    morpheus,
    plotly,
    RColorBrewer,
    scales,
    stringr,
    vegan
Suggests:
    testthat
Encoding: UTF-8
LazyData: true
Roxygen: list(old_usage=TRUE)
RoxygenNote: 7.1.1
```


##  R topics documented:
-   [datavis16s-package](#datavis16s-package)
-   [Exported](#exported)
    -   [adivboxplot](#adivboxplot)
    -   [allgraphs](#allgraphs)
    -   [morphheatmap](#morphheatmap)
    -   [pcoaplot](#pcoaplot)
    -   [rarefactioncurve](#rarefactioncurve)
    -   [readindata](#readindata)
    -   [trygraphwrapper](#trygraphwrapper)
-   [Internal](#internal)
    -   [amp_rarecurvefix](#amp_rarecurvefix)
    -   [filterlowabund](#filterlowabund)
    -   [gridCode](#gridcode)
    -   [highertax](#highertax)
    -   [log10scale](#log10scale)
    -   [logoutput](#logoutput)
    -   [plotlyGrid](#plotlygrid)
    -   [print_ampvis2](#print_ampvis2)
    -   [read_biom](#read_biom)
    -   [save_fillhtml](#save_fillhtml)
    -   [shortnames](#shortnames)
    -   [subsetamp](#subsetamp)

## datavis16s-package

datavis16s: A package for Nephele 16S pipeline visualization

## Exported

### adivboxplot

Alpha diversity boxplot

**Description**

Plots plotly boxplot of shannon diversity and Chao species richness. If
sampling depth is NULL, rarefies OTU table to the minimum readcount of
any sample. If this is low, then the plot will fail.

**Usage**

``` r
adivboxplot(datafile, outdir, mapfile, amp = NULL, sampdepth = NULL,
  colors = NULL, cats = NULL, filesuffix = NULL, ...)
```

**Arguments**

| Argument     | Description                                                                                                                                                  |
|--------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `datafile`   | full path to input OTU file                                                                                                                                  |
| `outdir`     | full path to output directory                                                                                                                                |
| `mapfile`    | full path to map file                                                                                                                                        |
| `amp`        | ampvis2 object. may be specified instead of mapfile and datafile                                                                                             |
| `sampdepth`  | sampling depth. see details.                                                                                                                                 |
| `colors`     | colors to use for plots                                                                                                                                      |
| `cats`       | categories/columns in mapping file to use as groups. If NULL (default), will use all columns starting with TreatmentGroup to (but not including) Description |
| `filesuffix` | (Optional) suffix for output filename                                                                                                                        |
| `...`        | other parameters to pass to [readindata](#readindata)                                                                                                        |

**Details**

If `sampdepth` is NULL, the sampling depth is set to the size of the
smallest sample.

**Value**

Save alpha diversity boxplots to outdir.

**Source**

[graphs.R](../R/graphs.R)

### allgraphs

Pipeline function

**Description**

Make all 4 types of graphs

**Usage**

``` r
allgraphs(datafile, outdir, mapfile, sampdepth = 10000, ...)
```

**Arguments**

| Argument    | Description                                                                             |
|-------------|-----------------------------------------------------------------------------------------|
| `datafile`  | full path to input OTU file (biom or txt file see [readindata](#readindata) for format) |
| `outdir`    | full path to output directory                                                           |
| `mapfile`   | full path to map file                                                                   |
| `sampdepth` | sampling depth. default: 10000                                                          |
| `...`       | other parameters to pass to [readindata](#readindata)                                   |

**Value**

graphs are saved to outdir. See [user doc](../doc/user_doc.md).

This value is used to remove samples before for alpha diversity and PCoA
plots. Also, to rarefy OTU table for the alpha diversity and Bray-Curtis
distance PCoA.

**Source**

[graphs.R](../R/graphs.R)

### morphheatmap

Morpheus heatmap

**Description**

Creates heatmaps using Morpheus R API
<https://software.broadinstitute.org/morpheus/> . The heatmaps are made
using relative abundances.

**Usage**

``` r
morphheatmap(datafile, outdir, mapfile, amp = NULL, sampdepth = NULL,
  rarefy = FALSE, filter_level = NULL, taxlevel = c("seq"),
  colors = NULL, rowAnnotations = NULL, force = FALSE,
  filesuffix = NULL, ...)
```

**Arguments**

| Argument         | Description                                                                                                                                        |
|------------------|----------------------------------------------------------------------------------------------------------------------------------------------------|
| `datafile`       | full path to input OTU file (biom or see [readindata](#readindata) )                                                                               |
| `outdir`         | full path to output directory                                                                                                                      |
| `mapfile`        | full path to mapping file                                                                                                                          |
| `amp`            | (Optional) ampvis2 object. may be specified instead of mapfile and datafile                                                                        |
| `sampdepth`      | sampling depth                                                                                                                                     |
| `rarefy`         | Logical. Rarefy the OTU table if sampdepth is specified.                                                                                           |
| `filter_level`   | minimum abundance to show in the heatmap                                                                                                           |
| `taxlevel`       | vector of taxonomic levels to graph. must be subset of c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "seq"). See Details. |
| `colors`         | (Optional) color vector - length equal to number of TreatmentGroups in mapfile                                                                     |
| `rowAnnotations` | (Optional) Row annotations to be used in addition to taxonomy.                                                                                     |
| `force`          | Force "seq" level heatmap to be made even if number of seqs is greater than 2000. ' See Details.                                                   |
| `filesuffix`     | (Optional) suffix for output filename                                                                                                              |
| `...`            | parameters to pass to [readindata](#readindata)                                                                                                    |

**Details**

For the `taxlevel` parameter, each level is made into a separate
heatmap. "seq" makes the heatmap with no collapsing of taxonomic levels
if there are fewer than 2000 ASVs/OTUs. Otherwise, Species level is made
instead.

**Value**

Saves heatmaps to outdir.

**Examples**

``` r
## Not run:
 morphheatmap(datafile="OTU_table.txt", outdir="outputs/graphs", mapfile="mapfile.txt",
sampdepth = 25000, taxlevel = c("Family", "seq"), tsvfile=TRUE) 
## End(Not run)
```

**Source**

[graphs.R](../R/graphs.R)

### pcoaplot

PCoA plots

**Usage**

``` r
pcoaplot(datafile, outdir, mapfile, amp = NULL, sampdepth = NULL,
  distm = "binomial", filter_species = 0.1, rarefy = FALSE,
  colors = NULL, filesuffix = NULL, ...)
```

**Arguments**

| Argument         | Description                                                                                                                                                                      |
|------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `datafile`       | full path to input OTU file (biom or see [readindata](#readindata) )                                                                                                             |
| `outdir`         | full path to output directory                                                                                                                                                    |
| `mapfile`        | full path to map file                                                                                                                                                            |
| `amp`            | ampvis2 object. may be specified instead of mapfile and datafile                                                                                                                 |
| `sampdepth`      | sampling depth                                                                                                                                                                   |
| `distm`          | distance measure for PCoA. any that are supported by [amp_ordinate](https://madsalbertsen.github.io/ampvis2/reference/amp_ordinate.html) except for unifrac, wunifrac, and none. |
| `filter_species` | Remove low abundant OTU's across all samples below this threshold in percent. Setting this to 0 may drastically increase computation time.                                       |
| `rarefy`         | Logical. Rarefy the OTU table if sampdepth is specified.                                                                                                                         |
| `colors`         | (Optional) color vector - length equal to number of TreatmentGroups in mapfile                                                                                                   |
| `filesuffix`     | (Optional) suffix for output filename                                                                                                                                            |
| `...`            | parameters to pass to [readindata](#readindata)                                                                                                                                  |

**Value**

Saves pcoa plots to outdir.

**Source**

[graphs.R](../R/graphs.R)

### rarefactioncurve

Make rarefaction curve graph

**Usage**

``` r
rarefactioncurve(datafile, outdir, mapfile, amp = NULL, colors = NULL,
  cat = "TreatmentGroup", stepsize = 1000, ...)
```

**Arguments**

| Argument   | Description                                                                                         |
|------------|-----------------------------------------------------------------------------------------------------|
| `datafile` | full path to input OTU file (biom or see [readindata](#readindata) )                                |
| `outdir`   | full path to output directory                                                                       |
| `mapfile`  | full path mapping file                                                                              |
| `amp`      | (Optional) ampvis2 object. may be specified instead of mapfile and datafile                         |
| `colors`   | (Optional) color vector - length equal to number of TreatmentGroups in mapfile                      |
| `cat`      | Category/column in mapping file by which to color the curves in the graph. (default TreatmentGroup) |
| `stepsize` | for rarefaction plotting.                                                                           |
| `...`      | parameters to pass to [readindata](#readindata)                                                     |

**Value**

Saves rarefaction curve plot to output directory.

**Source**

[graphs.R](../R/graphs.R)

### readindata

Read in data

**Usage**

``` r
readindata(datafile, mapfile, tsvfile = FALSE, mincount = 10)
```

**Arguments**

| Argument   | Description                                                                                     |
|------------|-------------------------------------------------------------------------------------------------|
| `datafile` | full path to input data file. must be either biom file or tab delimited text file. See details. |
| `mapfile`  | full path to mapfile. must contain SampleID, TreatmentGroup, and Description columns            |
| `tsvfile`  | Logical. Is datafile a tab-delimited text file? See details.                                    |
| `mincount` | minimum number of reads                                                                         |

**Details**

datafile may be either biom file or text file. If text file, it should
have ampvis2 OTU table format
<https://madsalbertsen.github.io/ampvis2/reference/amp_load.html#the-otu-table>.
If the number of reads is less than mincount, the function will give an
error, as we cannot make graphs with so few counts.

**Value**

ampvis2 object

**Source**

[graphs.R](../R/graphs.R)

### trygraphwrapper

Wrapper for any graph function

**Description**

This is a wrapper for any of the graph functions meant to be called
using rpy2 in python.

**Usage**

``` r
trygraphwrapper(datafile, outdir, mapfile, FUN, logfilename = "logfile.txt",
  info = TRUE, tsvfile = FALSE, ...)
```

**Arguments**

| Argument      | Description                                                                                           |
|---------------|-------------------------------------------------------------------------------------------------------|
| `datafile`    | full path to input OTU file (biom or txt, see [readindata](#readindata) for format of txt file)       |
| `outdir`      | output directory for graphs                                                                           |
| `mapfile`     | full path to map file                                                                                 |
| `FUN`         | character string. name of function you would like to run. can be actual function object if run from R |
| `logfilename` | logfilename                                                                                           |
| `info`        | print sessionInfo to logfile                                                                          |
| `tsvfile`     | Is datafile a tab-delimited text file? Default FALSE                                                  |
| `...`         | parameters needed to pass to FUN                                                                      |

**Value**

Returns 0 if FUN succeeds and stops on error. In rpy2, it will throw
rpy2.rinterface.RRuntimeError.

**Examples**

``` r
## Not run:
 # example with no optional arguments for running allgraphs
trygraphwrapper("/path/to/outputs/out.biom", "/path/to/outputs/",
"/path/to/inputs/mapfile.txt", 'allgraphs')

# example with sampdepth argument for running allgraphs
trygraphwrapper("/path/to/outputs/out.biom", "/path/to/outputs/",
"/path/to/inputs/mapfile.txt", 'allgraphs', sampdepth=30000)


# example with optional argument sampdepth and tsv file
trygraphwrapper("/path/to/outputs/OTU_table.txt", "/path/to/outputs/",
"/path/to/inputs/mapfile.txt", 'allgraphs', sampdepth = 30000, tsvfile=TRUE)

# example of making heatmap with optional arguments
trygraphwrapper("/path/to/outputs/taxa_species.biom", "/path/to/outputs",
"/path/to/inputs/mapfile.txt", 'morphheatmap', sampdepth = 30000, filter_level=0.01,
taxlevel=c("Family", "seq")) 
## End(Not run)
```

**Source**

[graphs.R](../R/graphs.R)

## Internal

### amp_rarecurvefix

Rarefaction curve

**Description**

This function replaces the ampvis2 function amp_rarecurve to fix
subsampling labeling bug in vegan

**Usage**

``` r
amp_rarecurvefix(data, stepsize = 1000, color_by = NULL)
```

**Arguments**

| Argument   | Description                                                                                  |
|------------|----------------------------------------------------------------------------------------------|
| `data`     | (required) Data list as loaded with amp_load.                                                |
| `stepsize` | Step size for the curves. Lower is prettier but takes more time to generate. (default: 1000) |
| `color_by` | Color curves by a variable in the metadata.                                                  |

**Value**

A ggplot2 object.

**Source**

[utilities.R](../R/utilities.R)

### filterlowabund

Filter low abundant taxa

**Usage**

``` r
filterlowabund(amp, level = 0.01, persamp = 0, abs = FALSE, toptaxa = NULL)
```

**Arguments**

| Argument  | Description                                                                                       |
|-----------|---------------------------------------------------------------------------------------------------|
| `amp`     | ampvis2 object                                                                                    |
| `level`   | level at which to filter                                                                          |
| `persamp` | percent of samples which must have taxa in common                                                 |
| `abs`     | is level an absolute count? if false, will use level as relative percent.                         |
| `toptaxa` | number of seqvar to include sorted by max count across all samples; if NULL all will be included. |

**Value**

filtered ampvis2 object

**Source**

[utilities.R](../R/utilities.R)

### gridCode

Format plotly grid code

**Description**

Format data according to here: <https://plot.ly/export/>

**Usage**

``` r
gridCode(data)
```

**Arguments**

| Argument | Description                  |
|----------|------------------------------|
| `data`   | data to populate plotly grid |

**Value**

list of 2 values:

-   `html` html for plotly export link  
-   `javascript` js function for exporting data

**Source**

[plotlyGrid.R](../R/plotlyGrid.R)

### highertax

return tables at higher tax level

**Usage**

``` r
highertax(amp, taxlevel)
```

**Arguments**

| Argument   | Description                                   |
|------------|-----------------------------------------------|
| `amp`      | ampvis2 object                                |
| `taxlevel` | taxonomic level at which to sum up the counts |

**Value**

ampvis2 object with otu table and taxa summed up to the taxlevel

**Source**

[utilities.R](../R/utilities.R)

### log10scale

Log base 10 + 1 scale

**Description**

Transformation which computes `log10(x+1)` scale

**Usage**

``` r
log10p_trans()
```

**Details**

`log10p` is for use with ggplot2 `trans` argument in scale function.

**Value**

`log10p` returns a scales tranformation object

### logoutput

write log output

**Description**

Prints time along with log message.

**Usage**

``` r
logoutput(c, bline = 0, aline = 0, type = NULL)
```

**Arguments**

| Argument | Description                                           |
|----------|-------------------------------------------------------|
| `c`      | String. Log message/command to print.                 |
| `bline`  | Number of blank lines to precede output.              |
| `aline`  | Number of blank lines to follow output.               |
| `type`   | String. Must be one of "WARNING", or "ERROR" or NULL. |

**Source**

[utilities.R](../R/utilities.R)

### plotlyGrid

Add Plotly data export to Plotly graph

**Description**

All functions create an output html plot with link which sends the data
to a grid in the plotly chart studio.

`plotlyGrid` takes in a ggplot or plotly object and creates an output
html plotly plot.

**Usage**

``` r
plotlyGrid(pplot, filename, data = NULL, title = NULL, outlib = "lib")
```

**Arguments**

| Argument   | Description                                                                                                            |
|------------|------------------------------------------------------------------------------------------------------------------------|
| `pplot`    | plotly or ggplot object                                                                                                |
| `filename` | output filename (fullpath)                                                                                             |
| `data`     | data frame to export to plotly grid (optional for plotlyGrid)                                                          |
| `title`    | title of html page                                                                                                     |
| `outlib`   | (Optional) name of external lib directory for non-selfcontained html. Useful for multiple graphs sharing the same lib. |

**Value**

html plot is saved to filename. external libraries are saved to outlib
in same directory as filename. Invisibly returns the plotly html widget.

**Source**

[plotlyGrid.R](../R/plotlyGrid.R)

### print_ampvis2

Print ampvis2 object summary

**Description**

This is a copy of the internal ampvis2 function print.ampvis2. CRAN does
not allow ':::' internal calling of function in package.

**Usage**

``` r
print_ampvis2(data)
```

**Arguments**

| Argument | Description    |
|----------|----------------|
| `data`   | ampvis2 object |

**Value**

Prints summary stats about ampvis2 object

**Source**

[utilities.R](../R/utilities.R)

### read_biom

biomformat read_biom

**Description**

This function replaces the biomformat function read_biom to deal with
reading in crappy hdf5 biom file.

**Usage**

``` r
read_biom(biom_file)
```

**Arguments**

| Argument    | Description          |
|-------------|----------------------|
| `biom_file` | input biom file name |

**Value**

biom object

### save_fillhtml

Save an HTML object to a file

**Usage**

``` r
save_fillhtml(html, file, background = "white", libdir = "lib", bodystyle = "")
```

**Arguments**

| Argument     | Description                       |
|--------------|-----------------------------------|
| `html`       | HTML content to print             |
| `file`       | File to write content to          |
| `background` | Background color for web page     |
| `libdir`     | Directory to copy dependencies to |
| `bodystyle`  | html style string                 |

**Value**

save html to file

**Source**

[plotlyGrid.R](../R/plotlyGrid.R)

### shortnames

shortnames for taxonomy

**Usage**

``` r
shortnames(taxtable)
```

**Arguments**

| Argument   | Description                                       |
|------------|---------------------------------------------------|
| `taxtable` | taxonomy table object from ampvis2 object amp$tax |

**Value**

data.frame taxonomy table object like ampvis2 amp$tax. taxonomy names
are sanitized and formatted to be a bit nicer.

**Source**

[utilities.R](../R/utilities.R)

### subsetamp

Subset and rarefy OTU table.

**Description**

Subset and/or rarefy OTU table.

**Usage**

``` r
subsetamp(amp, sampdepth = NULL, rarefy = FALSE, printsummary = T,
  outdir = NULL, ...)
```

**Arguments**

| Argument       | Description                                                                                                                          |
|----------------|--------------------------------------------------------------------------------------------------------------------------------------|
| `amp`          | ampvis2 object                                                                                                                       |
| `sampdepth`    | sampling depth. See details.                                                                                                         |
| `rarefy`       | rarefy the OTU table in addition to subsetting                                                                                       |
| `printsummary` | Logical. print ampvis2 summary of OTU table                                                                                          |
| `outdir`       | Output directory. If not null, and samples are removed from amp, the sample names will be output to outdir/samples_being_ignored.txt |
| `...`          | other parameters to pass to amp_subset_samples                                                                                       |

**Details**

`sampdepth` will be used to filter out samples with fewer than this
number of reads. If rarefy is TRUE, then it will also be used as the
depth at which to subsample using vegan function rrarefy.

**Value**

ampvis2 object

**Source**

[graphs.R](../R/graphs.R)
