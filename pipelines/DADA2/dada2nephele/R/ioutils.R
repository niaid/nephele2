
#' ASV sequence names
#'
#' @description DADA2 uses the actual ASV sequences as names for sequence and
#' taxonomy tables which is unwieldy for the output.  We make simpler new names
#' seq1, seq2, ...
#'
#' \emph{make_seq_names} makes the new names
#' @name make_seq_names
#'
#' @param seqtab OTU/DADA2 Sequence table
#' @param type one of: "simple" renaming seq1, seq2, ...
#' or "md5" for hash like in QIIME 2
#'
#' @return \emph{make_seq_names} returns named vector where elements are names and the names are the
#' sequences
#'
#' @source [ioutils.R](../R/ioutils.R#L18)
#'
#' @importFrom digest digest
#'
make_seq_names <- function(seqtab, type="simple") {
  if (type == "md5") {
    seqnames = sapply(colnames(seqtab), digest, algo="md5")
  } else if (type == "simple") {
    seqnames = paste0("seq", 1:ncol(seqtab))
  } else {
    stop("type must be one of {'simple', 'md5'}")
  }
  names(seqnames) = colnames(seqtab)
  return(seqnames)
}

#' @description \emph{replace_names} takes in old names and returns new ones.
#' @rdname make_seq_names
#'
#' @param tabnames character of original names
#' @param seq_names named vector with elements being new names and names from
#'  \code{tabnames}, e.g. output of \emph{make_seq_names}
#'
#' @return \emph{replace_names} returns (unnamed) vector of new names
#'
replace_names <- function(tabnames, seq_names) {
  newnames = unname(seq_names[tabnames])
  if (any(is.na(newnames))) stop("some sequences are not found in seq_names")
  return(newnames)
}


#' dada2output
#'
#' @description Utilities to convert dada2 sequence tables to output files.
#' @name dada2output
#'
#' @param otu OTU table.
#' @param tax Taxonomy table
#' @param metadata Metadata table (Optional).
#' @param filename Output filename
#'
#' @source [ioutils.R](../R/ioutils.R#L60)
#'
NULL

#' @description dada2biom makes valid biom object.
#'
#' @rdname dada2output
#'
#' @return dada2biom returns a biom object.
#'
#'
dada2biom <- function(otu, tax, metadata=NULL) {
  tax[is.na(tax)] <- "none"
  return(make_biom(t(otu), sample_metadata = metadata, observation_metadata = tax))
}


#' @rdname dada2output
#'
#' @description dada2text writes tab separated OTU table to filename.
#'
dada2text <- function(otu, tax, filename) {
  otutable <- cbind(t(otu), tax)
  write.table(otutable, file = filename, sep = "\t", quote = FALSE, na = "", col.names = NA)
}

#' @rdname dada2output
#'
#' @description dada2output writes tab separated taxonomy file in format
#' suitable for importing into qiime2
#'
dada2taxonomy <- function(tax, filename) {
  tax[is.na(tax)] <- "none"
  newtax <- data.frame(`Feature ID` = row.names(tax), check.names = FALSE)
  newtax$Taxon <- apply(tax, 1, function(x) paste(x, collapse="; "))
  write.table(newtax, file=filename, sep = "\t", quote=FALSE, na="", row.names = F)
}
