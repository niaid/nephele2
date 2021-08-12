#!/usr/bin/env Rscript

#' [Broad morpheus R lib](https://github.com/cmap/morpheus.R)
library(morpheus)
library(htmlwidgets)

#' Take in 2 arguments:
#'
#' 1. heatmap input *path_abun_unstrat_descrip.tsv*
#' 2. mapping file
#'
#' args <- commandArgs(trailingOnly = T)
#' heatmapinput <- args[1]
#' mapping <- args[2]
heatmapinput <- snakemake@input[[1]]
mapping <- snakemake@input[[2]]

#' read in files
input <- read.delim(heatmapinput, check.names = F)
metadata <- read.delim(mapping, check.names = F)
colnames(metadata) <- gsub("#SampleID", "SampleID", colnames(metadata))


#' make matrix for heatmap, **mat**, and data frame for row annotations, **rowmet**
rowmet = input[,1:2]
mat = input[,3:ncol(input)]
row.names(mat) <- rowmet$pathway

#' replace 0s as missing values and provide custom distance function **mydist** that
#' transforms NA back to 0
v <- which(mat == 0, arr.ind = T)
mat[v] <- NA
mydist <- function(x) {
  w <- which(is.na(x), arr.ind = T, useNames = F)
  x[w] <- 0
  return(dist(x))
}
row.names(mat) <- rowmet$pathway

#' make sure SampleID is character, not numeric.
#' check if metadata matches input file, make data frame **colmet** for column annotations,
#' and make list **coldisplay** to display column annotations
metadata$SampleID <- as.character(metadata$SampleID)
samples <- colnames(mat)
coldisplay <- list()

if (all(samples %in% metadata$SampleID) &&
    length(setdiff(colnames(metadata),c("SampleID", "ForwardFastqFile", "ReverseFastqFile"))) > 0) {

  colmet <- metadata[ which(metadata$SampleID %in% samples),
             setdiff(colnames(metadata), c("ForwardFastqFile", "ReverseFastqFile")) ]
  mat <- mat[, colmet$SampleID]

  ## display list
  coldisplay <- lapply(colnames(colmet)[2:ncol(colmet)], function(x) {
      return(list(field = x,
                  display = list('color'),
                  highlightMatchingValues = TRUE))
    })

} else {
  colmet <- NULL
}

#' add id to column display
coldisplay[[length(coldisplay) + 1]] <- list(field='id', display=list('text'))

#' save plain heatmap
mhmap <- morpheus(mat, rowAnnotations = rowmet, columnAnnotations=colmet,
                  columns = coldisplay, dist=mydist)
htmlwidgets::saveWidget(mhmap, snakemake@output[[1]], selfcontained = F)


