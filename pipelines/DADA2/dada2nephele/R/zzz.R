#' dada2 pipeline for nephele2
#'
#' @description This package implements Nephele 2.x DADA2 pipeline.
#'
#' @docType package
#' @name dada2nephele
#' @aliases dada2nephele-package
#'
#' @import dada2
#'
NULL

#' global constants
#'
#' @rdname onLoad
#' @description Set global constants to programmatically write documentation.  See [user doc](user_doc.md)
#' for more complete descriptions.  See examples for the actual values.
#'
#' @note The dparams option is set to be a list of the individual values as follows:
#' \describe{
#'   \item{maxEE}{parameter for dada2::filterAndTrim}
#'   \item{truncQ}{parameter for dada2::filterAndTrim}
#'   \item{trimLeft}{parameter for dada2::filterAndTrim}
#'   \item{nbases}{parameter for dada2::learnErrors}
#'   \item{minOverlap}{parameter for dada2::mergePairs}
#'   \item{maxMismatch}{parameter for dada2::mergePairs}
#'   \item{trimOverhang}{parameter for dada2::mergePairs}
#'   \item{justConcatenate}{parameter for dada2::mergePairs}
#'   \item{outputfasta}{filename for sequence variant FASTA file}
#'   \item{minBoot}{parameter for dada2::assignTaxonomy}
#'   \item{biomfile}{filename for biom file based on output of dada2::assignTaxonomy}
#'   \item{speciesbiomfile}{filename for biom file based on output of dada2::addSpecies}
#'   \item{otutable}{filename for tab delimited otu table based on output of dada2::addSpecies}
#'   \item{biomsummary}{filename for text file containing summary of OTU table}
#'   \item{refdb}{database filename}
#'   \item{refdb_species}{species database filename}
#'   \item{min_seq_length}{minimum length of denoised sequences to be used for taxonomic assignment.}
#'   \item{taxmethod}{method of taxonomic assignment}
#'   \item{taxtable}{filename for taxonomy table output}
#' }
#'
#' @examples
#' getOption("dparams")
#'
#' @source [onLoad:zzz.R](../R/zzz.R#L45)
#'
.onLoad <- function(libname, pkgname) {
  ## left trim parameter
  # trims <- data.frame(region=c("V1V3", "V4"), forward=c(20,19), reverse=c(17,20))
  # options(trims=trims)
  dparams <- list(nbases=1e8,
                  maxEE=5,
                  truncQ=4,
                  truncLen=0,
                  minOverlap=12,
                  maxMismatch=0,
                  justConcatenate=FALSE,
                  minBoot=80,
                  trimLeft=0,
                  trimOverhang=FALSE,
                  outputfasta='seq.fasta',
                  biomfile='taxa.biom',
                  otutable='OTU_table.txt',
                  biomsummary="otu_summary_table.txt",
                  refdb='dada2_silva_v132/silva_nr_v132_train_set.fa',
                  refdb_species='dada2_silva_v132/silva_species_assignment_v132.fa',
                  min_seq_length=75,
                  taxmethod='rdp',
                  taxtable='taxonomy_table.txt')
  options(dparams=dparams)
}
