################################################################################
#' Create a \link[biomformat]{biom-class}
#' from \code{\link{matrix-class}}
#' or \code{\link{data.frame}}.
#'
#' @note This code is forked from version 1.5 of the \href{https://github.com/joey711/biomformat-oldfork}{biomformat R library}.
#'
#' @description This function creates a valid instance of the \code{\link[biomformat]{biom-class}}
#' from standard base-R objects like
#' \code{\link{matrix-class}} or \code{\link{data.frame}}.
#' The object returned by this function is appropriate for writing to
#' a \code{.biom} file using the \code{\link{write_biom}} function.
#' The sparse biom-format is not (yet) supported.
#'
#' @param data (Required).
#'  \code{\link{matrix-class}} or \code{\link{data.frame}}.
#'  A contingency table.
#'  Observations / features / OTUs / species are rows,
#'  samples / sites / libraries are columns.
#'
#' @param sample_metadata (Optional).
#'  A \code{\link{matrix-class}} or \code{\link{data.frame}}
#'  with the number of rows equal to the number of samples in \code{data}.
#'  Sample covariates associated with the count data.
#'  This should look like the table returned by
#'  \code{\link[biomformat]{sample_metadata}} on a valid instance
#'  of the \code{\link[biomformat]{biom-class}}.
#'
#' @param observation_metadata (Optional).
#'  A \code{\link{matrix-class}} or \code{\link{data.frame}}
#'  with the number of rows equal to the number of
#'  features / species / OTUs / genes in \code{data}.
#'  This should look like the table returned by
#'  \code{\link[biomformat]{observation_metadata}} on a valid instance
#'  of the \code{\link[biomformat]{biom-class}}.
#'
#' @param id (Optional). Character string. Identifier for the project.
#'
#' @param matrix_element_type (Optional). Character string. Either 'int' or 'float'
#'
#' @param qiime_format (Optional). Logical. biom-format requires that observation
#' metadata be key, value pairs (or group,
#'  dataset for hd5).  For QIIME, there is only one pair with key be set to "taxonomy,"
#'  and the value must be the entire taxonomy table.  If FALSE, each column of observation
#'  metadata will be a separate key (to be consistent with sample metadata).
#'
#' @return An object of \code{\link[biomformat]{biom-class}}.
#'
#' @importFrom biomformat biom
#' @importFrom utils packageVersion
#' @importFrom methods as
#'
#' @source [make_biom in biomformat.R](../R/biomformat.R#L53)
#'
make_biom <- function(data, sample_metadata=NULL, observation_metadata=NULL, id=NULL, matrix_element_type="int", qiime_format=TRUE){
  # The observations / features / OTUs / rows "meta" data table
  if(!is.null(observation_metadata)){
    rows = mapply(list, SIMPLIFY=FALSE, id=as.list(rownames(data)),
                  metadata=lapply(seq_len(nrow(observation_metadata)), function(i) observation_metadata[i,]))
  } else {
    rows = mapply(list, id=as.list(rownames(data)), metadata=NA, SIMPLIFY=FALSE)
  }
  # The samples / sites / columns "meta" data table
  if(!is.null(sample_metadata)){
    columns = mapply(list, SIMPLIFY=FALSE, id=as.list(colnames(data)),
                     metadata=lapply(seq_len(nrow(sample_metadata)), function(i) sample_metadata[i,]))
  } else {
    columns = mapply(list, id=as.list(colnames(data)), metadata=NA, SIMPLIFY=FALSE)
  }
  # Convert the contingency table to a list
  datalist = as.list(as.data.frame(as(t(data), "matrix")))
  names(datalist) <- NULL
  # Define the list, instantiate as biom-format, and return
  # (Might eventually expose some of these list elements as function arguments)
  format_url = "http://biom-format.org/documentation/format_versions/biom-1.0.html"
  biomout <- biom(list(id=id,
                       format = "Biological Observation Matrix 1.0.0",
                       format_url = format_url,
                       type = "OTU table",
                       generated_by = sprintf("biomformat %s", packageVersion("biomformat")),
                       date = strftime(Sys.time(), format="%Y-%m-%dT%H:%M:%S"),
                       matrix_type = "dense",
                       matrix_element_type = matrix_element_type,
                       shape = dim(data),
                       rows = rows,
                       columns = columns,
                       data = datalist))
  if (!is.null(observation_metadata)) {
    if (qiime_format) {
      biomout$rows <- lapply(biomout$rows, function(x) { mm <- x$metadata; x$metadata <- NULL; x$metadata$taxonomy <- mm; x })
    } else {
      biomout$rows <- lapply(biomout$rows, function(x) { x$metadata <- as.list(x$metadata); return(x) })
    }
  }

  if (!is.null(sample_metadata)) {
    biomout$columns <- lapply(biomout$columns, function(x) { x$metadata <- as.list(x$metadata); return(x) })
  }

  return(biomout)
}

################################################################################
#' Write a biom-format v1 file, returning a \link[biomformat]{biom-class}.
#'
#' @note This code is forked from version 1.5 of the \href{https://github.com/joey711/biomformat-oldfork}{biomformat R library}.
#'
#' @param biom (Required). A biom object.
#'
#' @param biom_file (Required). A character string indicating the
#'  file location of the biom formatted file. This is a JSON formatted file
#'  specific to biological datasets.
#'  The format is formally defined at
#'  \href{http://biom-format.org/documentation/biom_format.html}{the biom-format definition}
#'
#' @param pretty logical; Should biom output be pretty printed?
#'
#' @return Nothing. The first argument, \code{x}, is written to a file.
#'
#'
#' @importFrom jsonlite toJSON write_json
#'
#' @source [write_biom in biomformat.R](../R/biomformat.R#L123)
#'
write_biom <- function(biom, biom_file, pretty=FALSE) {
  ## deal with case of only one sample
  biom$data <- lapply(biom$data, I)
  write_json(biom, path=biom_file, pretty=pretty, auto_unbox = TRUE, na = "null")
}
