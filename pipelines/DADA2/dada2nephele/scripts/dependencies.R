#!/usr/bin/env Rscript

## Nephele is using R 3.5.x. The following script only works with this version!
if (as.numeric(version$major) != 3 || as.numeric(version$minor) >= 6 || as.numeric(version$minor) < 5) stop("This package only works with R version 3.5.x.  The version you are using is not compatible.")

options(repos=c(CRAN="https://cran.rstudio.com/"))

## Bioconductor 3.8
if (!("BiocManager" %in% installed.packages()) || packageVersion("BiocManager") != "1.30.4") install.packages("BiocManager")

## bioconductor packages
if (!"biomformat" %in% installed.packages() || packageVersion("biomformat") != "1.10.1" ) BiocManager::install('biomformat', version="3.8")
if (!"dada2" %in% installed.packages() || packageVersion("dada2") != "1.10.1" ) BiocManager::install('dada2', version="3.8")
if (!"DECIPHER" %in% installed.packages() || packageVersion("DECIPHER") != "2.10.1" ) BiocManager::install('DECIPHER', version="3.8")

## devtools
if (!"devtools" %in% installed.packages()) install.packages('devtools')

## dada2nephele
## May have to change path to dada2nephele
devtools::install_local("/usr/local/src/nephele2/pipelines/DADA2/dada2nephele", dependencies=TRUE, force=TRUE)
