#!/usr/bin/env Rscript

## See packer install for how we built AMI using conda: https://github.niaid.nih.gov/bcbb/nephele2-amis/blob/master/dada2_qiime2_2020_11/install.sh
## Below is how to mostly replicate in R

## Nephele is using R 4.0.x. The following script only works with this version!
if (as.numeric(version$major) != 4 || as.numeric(version$minor) >= 0 ) stop("This script only works with R version 4.0.x.  The version you are using is not compatible.  This package has only been tested with R 4.0.3 and Bioconductor 3.12.")

options(repos=c(CRAN="https://cran.rstudio.com/"))

## Bioconductor 3.12
if (!("BiocManager" %in% installed.packages())) install.packages("BiocManager")

## bioconductor packages
if (!"biomformat" %in% installed.packages() || packageVersion("biomformat") != "1.18" ) BiocManager::install('biomformat', version="3.12")
if (!"dada2" %in% installed.packages() || packageVersion("dada2") != "1.18" ) BiocManager::install('dada2', version="3.12")
if (!"DECIPHER" %in% installed.packages() || packageVersion("DECIPHER") != "2.18" ) BiocManager::install('DECIPHER', version="3.12")

## devtools
if (!"remotes" %in% installed.packages()) install.packages('remotes')

## dada2nephele
## May have to change path to dada2nephele
remotes::install_local("/usr/local/src/nephele2/pipelines/DADA2/dada2nephele",
                       dependencies=TRUE, force=TRUE)
