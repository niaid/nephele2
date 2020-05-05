#!/usr/bin/env Rscript

## Help for functions here:  https://githuburl.whatever/nephele2/blob/master/pipelines/DADA2/dada2nephele/doc/Reference_Manual_dada2nephele.md
## Main function dada2compute.


## Load dada2nephele package
library(dada2nephele)

## Option dparams holds all default parameters for dada2 pipeline.  Some of them can be set by the user in Nephele, and these parameters are passed directly to dada2compute.
## Non-user parameters can be modified directly with dparams.
## Example
dparams <- getOption("dparams")
print(dparams)
dparams$truncQ <- c(2,2)
options(dparams=dparams)

## Main function is dada2compute.  Can do ?dada2compute or see function help (link above).
## Intermediate files get saved to the working directory.
## Final files get saved to outdir.
outdir <- "example/test"
setwd(outdir)
dada2compute("example/input", outdir , "example/testmapfile.txt", refdb="silva_nr_v128_train_set.fa.gz", refdb_species="silva_species_assignment_v128.fa.gz", trimLeft=c(19,20), trimOverhang=FALSE, nthread=30, chimera=FALSE)


