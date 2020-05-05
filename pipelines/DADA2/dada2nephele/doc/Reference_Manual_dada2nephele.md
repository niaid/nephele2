Package 'dada2nephele'
=======================


```
Type: Package
Package: dada2nephele
Title: dada2 package for nephele2
Version: 0.1.1
Authors@R (parsed):
    * OCICB NIAID/NIH <nephelesupport@nih.gov> [aut, cre]
Description: Implements the DADA2 pipeline for use in the Nephele 2.x
    environment with the python library rpy2.
URL:
    https://github.com/niaid/nephele2/tree/master/pipelines/DADA2/dada2nephele
Depends:
    methods
Imports:
    biomformat,
    Biostrings,
    dada2,
    DECIPHER,
    foreach,
    ggplot2,
    jsonlite,
    ShortRead,
    stringr
ByteCompile: true
Encoding: UTF-8
LazyData: true
Roxygen: list(old_usage = FALSE)
RoxygenNote: 7.1.0
```


##  R topics documented:
  - [Exported](#exported)
      - [dada2compute](#dada2compute)
  - [Internal](#internal)
      - [checktrimfiles](#checktrimfiles)
      - [comb](#comb)
      - [dada2nephele-package](#dada2nephele-package)
      - [dada2output](#dada2output)
      - [decipher\_assign](#decipher_assign)
      - [getN](#getn)
      - [logoutput](#logoutput)
      - [make\_biom](#make_biom)
      - [.onLoad](#onload)
      - [run\_cmd](#run_cmd)
      - [write\_biom](#write_biom)

## Exported

### dada2compute

run dada2 pipeline

**Description**

`dada2compute` will run the whole pipeline from within R. See
[.onLoad](#onload) examples for argument default values or
`getOption('dparams')` .

`trycomputewrapper` is a wrapper for `dada2compute` to be used with rpy2
to run the pipeline in python. It sets up the global R options and the
output to the logfile before running `dada2compute` .

**Usage**

``` r
dada2compute(
  datadir,
  outdir,
  mapfile,
  refdb = getOption("dparams")$refdb,
  refdb_species = getOption("dparams")$refdb_species,
  nthread = TRUE,
  chimera = FALSE,
  trimLeft = getOption("dparams")$trimLeft,
  trimOverhang = getOption("dparams")$trimOverhang,
  data_type = "PE",
  maxEE = getOption("dparams")$maxEE,
  truncQ = getOption("dparams")$truncQ,
  minBoot = getOption("dparams")$minBoot,
  truncLen = getOption("dparams")$truncLen,
  maxMismatch = getOption("dparams")$maxMismatch,
  justConcatenate = getOption("dparams")$justConcatenate,
  taxmethod = getOption("dparams")$taxmethod,
  band_size = NULL,
  homopolymer_gap_penalty = NULL
)
trycomputewrapper(datadir, outdir, mapfile, logfilename = "logfile.txt", ...)
```

**Arguments**

| Argument                  | Description                                                                                                                                                                                                                                                                                       |
| ------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `datadir`                 | Directory containing FASTQ files.                                                                                                                                                                                                                                                                 |
| `outdir`                  | Output directory.                                                                                                                                                                                                                                                                                 |
| `mapfile`                 | Mapping filename - needs to be tab separated and have, at minimum, SampleID, ForwardFastqFile, and ReverseFastqFile columns.                                                                                                                                                                      |
| `refdb`                   | Path to the reference database.                                                                                                                                                                                                                                                                   |
| `refdb_species`           | Path to the species reference database.                                                                                                                                                                                                                                                           |
| `nthread`                 | Number of processors for parallel steps.                                                                                                                                                                                                                                                          |
| `chimera`                 | Should chimeric sequences be removed? If primers are not                                                                                                                                                                                                                                          |
| `trimLeft`                | (Optional). In [dada2::filterAndTrim](#dada2::filterandtrim) , the number of nucleotides to remove from the start of each read, forward and reverse. The values should be chosen based on the lengths of primers used for sequencing. (Whole number vector of length 2, for forward and reverse). |
| `trimOverhang`            | (Optional). After merging paired end reads, trim sequence which overhangs the start of each read. If amplicons are shorter than read length, e.g. 16S V4 region, we suggest setting this to True. (Logical. Default: False).                                                                      |
| `data_type`               | 'SE' or 'PE'                                                                                                                                                                                                                                                                                      |
| `maxEE`                   | Integer. Remove reads with less than maxEE expected errors in [dada2::filterAndTrim](#dada2::filterandtrim)                                                                                                                                                                                       |
| `truncQ`                  | Integer. Truncate reads at first qual below truncQ in [dada2::filterAndTrim](#dada2::filterandtrim)                                                                                                                                                                                               |
| `minBoot`                 | Integer. Bootstrap value for [dada2::assignTaxonomy](#dada2::assigntaxonomy)                                                                                                                                                                                                                      |
| `truncLen`                | Integer. Truncation length for reads in [dada2::filterAndTrim](#dada2::filterandtrim) . Can be vector of length 1 or 2                                                                                                                                                                            |
| `maxMismatch`             | Integer. Maximum number of allowed mismatches in [dada2::mergePairs](#dada2::mergepairs)                                                                                                                                                                                                          |
| `justConcatenate`         | Logical. Should PE reads be concatenated instead of overlapped and merged.                                                                                                                                                                                                                        |
| `taxmethod`               | method for taxonomic assignment. either idtaxa or rdp.                                                                                                                                                                                                                                            |
| `band_size`               | Integer. Band size for [dada2::dada](#dada2::dada) , [dada2::setDadaOpt](#dada2::setdadaopt)                                                                                                                                                                                                      |
| `homopolymer_gap_penalty` | Integer. for [dada2::dada](#dada2::dada) & [dada2::setDadaOpt](#dada2::setdadaopt) . cost of gaps in homopolymer regions (\>=3 repeated bases)                                                                                                                                                    |
| `logfilename`             | Log file name (full path).                                                                                                                                                                                                                                                                        |
| `...`                     | parameters to pass to `dada2compute`                                                                                                                                                                                                                                                              |

**Value**

`trycomputewrapper` returns 0 if `dada2compute` succeeds and stops on
error. In rpy2, it will throw `rpy2.rinterface.RRuntimeError` .

**Source**

[dada2compute: computation.R](../R/computation.R#L366)

[trycomputewrapper: computation.R](../R/computation.R#L635)

## Internal

### checktrimfiles

Check files after filterAndTrim

**Usage**

``` r
checktrimfiles(A, filt.dir, trimlist)
```

**Arguments**

| Argument   | Description                                 |
| ---------- | ------------------------------------------- |
| `A`        | mapping data.frame                          |
| `filt.dir` | filtered data directory                     |
| `trimlist` | list of vectors, R1 and R2 of trimmed files |

**Value**

list of `A` , `trimr1` , and maybe `trimr2` with missing files/rows
removed. Stops on error if no trimmed files exist.

**Source**

[computation.R](../R/computation.R#L207)

### comb

Combine function foreach loop

**Description**

For a [foreach](#foreach) loop that returns a list of n items for each
iteration, this function combines all iterations into n different lists
- one for each item returned.

**Usage**

``` r
comb(x, ...)
```

**Arguments**

| Argument | Description                                             |
| -------- | ------------------------------------------------------- |
| `x`      | list of n items to append onto                          |
| `...`    | individual lists of n items to append onto n lists in x |

**Details**

Use this as the value of the `.combine` parameter in foreach. You will
need to also initiate the list of n lists for the `.init` foreach
parameter.

**Value**

list of lists named according to names of list from single iteration.

**Examples**

``` r
## Not run:
 `%do%` <- foreach::`%do%`
oper <- foreach::foreach(i=1:5, .combine='comb', .multicombine=TRUE,
.init=list(list(), list())) %do% {
list(i*2, i*3)
}

oper[[1]]
## End(Not run)
```

### dada2nephele-package

dada2nephele: dada2 package for nephele2

### dada2output

dada2output

**Description**

Utilities to convert dada2 sequence tables to output files.

dada2biom makes valid biom object.

dada2text writes tab separated OTU table to filename.

dada2output writes tab separated taxonomy file in format suitable for
importing into qiime2

**Usage**

``` r
dada2biom(otu, tax, metadata = NULL)
dada2text(otu, tax, filename)
dada2taxonomy(tax, filename)
```

**Arguments**

| Argument   | Description                |
| ---------- | -------------------------- |
| `otu`      | OTU table.                 |
| `tax`      | Taxonomy table             |
| `metadata` | Metadata table (Optional). |
| `filename` | Output filename            |

**Value**

dada2biom returns a biom object.

**Source**

[computation.R](../R/computation.R#L156)

### decipher\_assign

Assign taxonomy using DECIPHER

**Usage**

``` r
decipher_assign(refdb, seqtab, nthread)
```

**Arguments**

| Argument  | Description                    |
| --------- | ------------------------------ |
| `refdb`   | reference database .RData file |
| `seqtab`  | sequence table made by DADA2   |
| `nthread` | number of processors to use    |

**Value**

taxa table in DADA2 taxa format

### getN

Get number of reads

**Description**

get number of reads for each sample from a dada2 object (such as output
of [mergePairs](#mergepairs) ).

**Usage**

``` r
getN(x)
```

**Arguments**

| Argument | Description  |
| -------- | ------------ |
| `x`      | dada2 object |

**Value**

Integer vector

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
| -------- | ----------------------------------------------------- |
| `c`      | String. Log message/command to print.                 |
| `bline`  | Number of blank lines to precede output.              |
| `aline`  | Number of blank lines to follow output.               |
| `type`   | String. Must be one of "WARNING", or "ERROR" or NULL. |

**Source**

[computation.R](../R/computation.R#L83)

### make\_biom

Create a [biom-class](#biom-class) from [matrix-class](#matrix-class) or
[data.frame](#dataframe) .

**Description**

This function creates a valid instance of the [biom-class](#biom-class)
from standard base-R objects like [matrix-class](#matrix-class) or
[data.frame](#dataframe) . The object returned by this function is
appropriate for writing to a `.biom` file using the
[write\_biom](#writebiom) function. The sparse biom-format is not (yet)
supported.

**Usage**

``` r
make_biom(
  data,
  sample_metadata = NULL,
  observation_metadata = NULL,
  id = NULL,
  matrix_element_type = "int",
  qiime_format = TRUE
)
```

**Arguments**

| Argument               | Description                                                                                                                                                                                                                                                                                                                                          |
| ---------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `data`                 | (Required). [matrix-class](#matrix-class) or [data.frame](#dataframe) . A contingency table. Observations / features / OTUs / species are rows, samples / sites / libraries are columns.                                                                                                                                                             |
| `sample_metadata`      | (Optional). A [matrix-class](#matrix-class) or [data.frame](#dataframe) with the number of rows equal to the number of samples in `data` . Sample covariates associated with the count data. This should look like the table returned by [sample\_metadata](#samplemetadata) on a valid instance of the [biom-class](#biom-class) .                  |
| `observation_metadata` | (Optional). A [matrix-class](#matrix-class) or [data.frame](#dataframe) with the number of rows equal to the number of features / species / OTUs / genes in `data` . This should look like the table returned by [observation\_metadata](#observationmetadata) on a valid instance of the [biom-class](#biom-class) .                                |
| `id`                   | (Optional). Character string. Identifier for the project.                                                                                                                                                                                                                                                                                            |
| `matrix_element_type`  | (Optional). Character string. Either 'int' or 'float'                                                                                                                                                                                                                                                                                                |
| `qiime_format`         | (Optional). Logical. biom-format requires that observation metadata be key, value pairs (or group, dataset for hd5). For QIIME, there is only one pair with key be set to "taxonomy," and the value must be the entire taxonomy table. If FALSE, each column of observation metadata will be a separate key (to be consistent with sample metadata). |

**Value**

An object of [biom-class](#biom-class) .

**Note**

This code is forked from version 1.4 of the [biomformat R
library](https://github.com/joey711/biomformat-oldfork) released under
[GPL-2](https://www.r-project.org/Licenses/GPL-2) .

### .onLoad

global constants

**Description**

Set global constants to programmatically write documentation. See [user
doc](user_doc.md) for more complete descriptions. See examples for the
actual values.

**Usage**

``` r
.onLoad(libname, pkgname)
```

**Note**

The dparams option is set to be a list of the individual values as
follows:

  - `maxEE` parameter for dada2::filterAndTrim
  - `truncQ` parameter for dada2::filterAndTrim
  - `trimLeft` parameter for dada2::filterAndTrim
  - `nbases` parameter for dada2::learnErrors
  - `minOverlap` parameter for dada2::mergePairs
  - `maxMismatch` parameter for dada2::mergePairs
  - `trimOverhang` parameter for dada2::mergePairs
  - `justConcatenate` parameter for dada2::mergePairs
  - `outputfasta` filename for sequence variant FASTA file
  - `minBoot` parameter for dada2::assignTaxonomy
  - `biomfile` filename for biom file based on output of
    dada2::assignTaxonomy
  - `speciesbiomfile` filename for biom file based on output of
    dada2::addSpecies
  - `otutable` filename for tab delimited otu table based on output of
    dada2::addSpecies
  - `biomsummary` filename for text file containing summary of OTU
    table
  - `refdb` database filename
  - `refdb_species` species database filename
  - `min_seq_length` minimum length of denoised sequences to be used for
    taxonomic assignment.
  - `taxmethod` method of taxonomic assignment
  - `taxtable` filename for taxonomy table output

**Examples**

``` r
getOption("dparams")
```

    ## $nbases
    ## [1] 100000000
    ##
    ## $maxEE
    ## [1] 5
    ##
    ## $truncQ
    ## [1] 4
    ##
    ## $truncLen
    ## [1] 0
    ##
    ## $minOverlap
    ## [1] 12
    ##
    ## $maxMismatch
    ## [1] 0
    ##
    ## $justConcatenate
    ## [1] FALSE
    ##
    ## $minBoot
    ## [1] 80
    ##
    ## $trimLeft
    ## [1] 20
    ##
    ## $trimOverhang
    ## [1] FALSE
    ##
    ## $outputfasta
    ## [1] "seq.fasta"
    ##
    ## $biomfile
    ## [1] "taxa.biom"
    ##
    ## $speciesbiomfile
    ## [1] "taxa_species.biom"
    ##
    ## $otutable
    ## [1] "OTU_table.txt"
    ##
    ## $biomsummary
    ## [1] "otu_summary_table.txt"
    ##
    ## $refdb
    ## [1] "dada2_silva_v132/silva_nr_v132_train_set.fa"
    ##
    ## $refdb_species
    ## [1] "dada2_silva_v132/silva_species_assignment_v132.fa"
    ##
    ## $min_seq_length
    ## [1] 75
    ##
    ## $taxmethod
    ## [1] "rdp"
    ##
    ## $taxtable
    ## [1] "taxonomy_table.txt"

**Source**

[onLoad:computation.R](../R/computation.R#L43)

### run\_cmd

wrap command in tryCatch

**Description**

(optionally) log cmd to output and evaluate cmd in the parent
environment, catching and parsing errors.

**Usage**

``` r
run_cmd(cmd, step = NULL, bline = 0, aline = 0, w2e = NA, log = T)
```

**Arguments**

| Argument | Description                                                                        |
| -------- | ---------------------------------------------------------------------------------- |
| `cmd`    | Character. Command string.                                                         |
| `step`   | (Optional). Step name to pass to error. Default will use cmd.                      |
| `bline`  | Number of blank lines to precede log output; parameter for [logoutput](#logoutput) |
| `aline`  | Number of blank lines to follow log output; parameter for [logoutput](#logoutput)  |
| `w2e`    | warning message to escalate to error.                                              |
| `log`    | log command to file                                                                |

**Value**

cmd will be evaluated in the parent environment, so return values in cmd
will be there.

**Source**

[computation.R](../R/computation.R#L113)

### write\_biom

Write a biom-format v1 file, returning a [biom-class](#biom-class) .

**Usage**

``` r
write_biom(biom, biom_file, pretty = FALSE)
```

**Arguments**

| Argument    | Description                                                                                                                                                                                                                                                                  |
| ----------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `biom`      | (Required). A biom object.                                                                                                                                                                                                                                                   |
| `biom_file` | (Required). A character string indicating the file location of the biom formatted file. This is a JSON formatted file specific to biological datasets. The format is formally defined at [the biom-format definition](http://biom-format.org/documentation/biom_format.html) |
| `pretty`    | logical; Should biom output be pretty printed?                                                                                                                                                                                                                               |

**Value**

Nothing. The first argument, `x` , is written to a file.

**Note**

This code is forked from version 1.4 of the [biomformat R
library](https://github.com/joey711/biomformat-oldfork) released under
[GPL-2](https://www.r-project.org/Licenses/GPL-2) .
