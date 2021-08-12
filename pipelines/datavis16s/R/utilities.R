#' write log output
#'
#' @description Prints time along with log message.
#'
#' @param c String. Log message/command to print.
#' @param bline Number of blank lines to precede output.
#' @param aline Number of blank lines to follow output.
#' @param type  String.  Must be one of "WARNING", or "ERROR" or NULL.
#'
#' @source [utilities.R](../R/utilities.R)
#'
logoutput <- function(c, bline = 0, aline = 0, type=NULL) {

  if (!is.null(type)) {
    if (!(type %in% c('WARNING', 'ERROR'))) {
      stop('type parameter must be one of "WARNING" or "ERROR" or NULL.')
    }
    type <- paste0(" - ", type)
  }

  lo <- paste0("[", Sys.time(), type, "] ", c)
  writeLines(c(rep("", bline), lo, rep("", aline)))

}


#' shortnames for taxonomy
#'
#' @param taxtable taxonomy table object from ampvis2 object amp$tax
#'
#' @return data.frame taxonomy table object like ampvis2 amp$tax.  taxonomy names
#' are sanitized and formatted to be a bit nicer.
#'
#' @source [utilities.R](../R/utilities.R)
#'
shortnames <- function(taxtable) {
  ## greengenes
  taxtable <- apply(taxtable, 2, function(x) { v <- grep("_unclassified$|^Unassigned$", x); x[v] <- NA; x })
  taxtable <- apply(taxtable, 2, function(x) { gsub("^[kpcofgs]__", "", x) })

  ## SILVA97
  taxtable <- apply(taxtable, 2, function(x) { gsub("^D_\\d+__", "", x) })
  uncnames <-  c("Other", "uncultured", "uncultured bacterium", "Ambiguous_taxa", "uncultured organism", "uncultured rumen bacterium")

  taxtable[which(taxtable %in% c("", "none"), arr.ind = TRUE)] <- NA


  ## Unclassified
  vd <- which(is.na(taxtable[,"Kingdom"]))
  taxtable[vd,"Kingdom"] <- paste(row.names(taxtable)[vd], "unclassified")


  ## Output Genus species instead of just species
  vn <- which(apply(taxtable,1, function(x) !(is.na(x["Species"]) | x["Species"] %in% uncnames | grepl(paste0("^",stringr::str_replace_all(x['Genus'], "(\\W)", "\\\\\\1")), x["Species"]))))
  if (length(vn) > 0) {
    taxtable[vn,"Species"] <- paste(taxtable[vn,"Genus"],taxtable[vn,"Species"] )
  }


  ## Generate newnames
  snamecol <- function(taxa) {
    sname <- function(x) {
      ## SILVA
      y <- suppressWarnings(min(which(x %in% uncnames)))
      if (y != "Inf"){
        return(paste(x[y-1], colnames(taxtable)[y-1], x[y]))
      }

      y <- suppressWarnings(min(which(is.na(x))))
      if (y == "Inf") {
        return(paste(x[length(x)]))
      }

      if ( y == 2 ) {
        return(paste(x[1]))
      }


      return(paste(x[y-1], colnames(shorttable)[y-1]))

    }

    tt <- which(colnames(taxtable) == taxa)
    shorttable <- taxtable[,1:tt]
    sn <- apply(shorttable, 1, sname )
    return(sn)
  }


  taxalist <- intersect(c("Phylum", "Class", "Order", "Family", "Genus", "Species"), colnames(taxtable))
  tlist <- lapply(taxalist, snamecol)
  newname_taxtable <- do.call(cbind.data.frame, tlist)
  newname_taxtable <- cbind.data.frame(taxtable[,"Kingdom"], newname_taxtable, taxtable[,"OTU"])
  colnames(  newname_taxtable ) <- colnames(taxtable)
  return(newname_taxtable)

}

#' return tables at higher tax level
#'
#' @param amp  ampvis2 object
#' @param taxlevel  taxonomic level at which to sum up the counts
#'
#' @return  ampvis2 object with otu table and taxa summed up to the taxlevel
#'
#' @importFrom data.table data.table .SD
#'
#' @source [utilities.R](../R/utilities.R)
#'
highertax <- function(amp, taxlevel) {

  otu <- as.matrix(amp$abund)
  tax <- amp$tax

  tc <- which(colnames(tax) == taxlevel)
  sn <- shortnames(tax)
  tax <- tax[,1:tc]
  tax[,tc] <- sn[,tc]

  taxcols <- colnames(tax)
  otucols <- colnames(otu)
  totu <- cbind.data.frame(otu, tax)
  dt <- data.table(totu)
  dt <- dt[, lapply(.SD, sum) , by = c(taxcols)]
  otu <- data.frame(dt[,otucols, with=FALSE], check.names = FALSE)
  tax <- data.frame(dt[,taxcols, with=FALSE])
  dupes <- tax[duplicated(tax[,tc]),tc]
  dupes <- which(tax[,tc] %in% dupes)
  tax[dupes,tc] <- paste(tax[dupes,tc-1], tax[dupes,tc])

  row.names(otu) <- tax[,tc]
  row.names(tax) <- tax[,tc]
  amp$abund <- as.data.frame(otu)
  amp$tax <- tax
  return(amp)
}

#' Filter low abundant taxa
#'
#' @param amp  ampvis2 object
#' @param level  level at which to filter
#' @param persamp  percent of samples which must have taxa in common
#' @param abs  is level an absolute count? if false, will use level as relative percent.
#' @param toptaxa number of seqvar to include sorted by max count across all samples;
#' if NULL all will be included.
#'
#' @return filtered ampvis2 object
#'
#' @source [utilities.R](../R/utilities.R)
#'
filterlowabund <- function(amp, level=0.01, persamp=0, abs=FALSE, toptaxa=NULL) {

  otu <- as.matrix(amp$abund)
  tax <- amp$tax
  if (abs) {
    mat <- otu
  } else {
    mat <- apply(otu, 2, function(x) { y <- sum(x); 100*x/y  } )
  }
  v <- which(mat < level & mat > 0, arr.ind = T)
  otu[v] <- 0
  w <- which(rowSums(otu) == 0)
  ps <- which(apply(otu, 1, function(x) length(which(x > 0))/ncol(otu)) < persamp/100)
  w <- unique(c(w, ps))
  if (length(w) > 0) {
    otu <- otu[-w,]
    tax <- tax[-w,]
  }

  if (!is.null(toptaxa)) {
    mv <- apply(otu, 1, max)
    mw <- order(mv, decreasing = T)[1:toptaxa]
    otu <- otu[mw,]
    tax <- tax[mw,]
  }

  amp$abund <- as.data.frame(otu)
  amp$tax <- tax
  return(amp)
}


#' Print ampvis2 object summary
#'
#' @param data ampvis2 object
#'
#' @return Prints summary stats about ampvis2 object
#'
#' @source [utilities.R](../R/utilities.R)
#'
#' @description  This is a copy of the internal ampvis2 function print.ampvis2.  CRAN does not allow
#' ':::' internal calling of function in package.
#'
#' @importFrom stats median
#'
print_ampvis2 <- function(data) {

  cat(class(data), "object with", length(data), "elements.\nSummary of OTU table:\n")
  print.table(c(
    Samples = as.character(ncol(data$abund)), OTUs = as.character(nrow(data$abund)),
    `Total#Reads` = as.character(sum(data$abund)), `Min#Reads` = as.character(min(colSums(data$abund))),
    `Max#Reads` = as.character(max(colSums(data$abund))),
    `Median#Reads` = as.character(median(colSums(data$abund))),
    `Avg#Reads` = as.character(round(mean(colSums(data$abund)),
      digits = 2
    ))
  ), justify = "right")
  cat("\nAssigned taxonomy:\n")
  print.table(c(Kingdom = paste(sum(nchar(data$tax$Kingdom) >
    3), "(", round(sum(nchar(data$tax$Kingdom) > 3) / nrow(data$abund),
    digits = 2
  ) * 100, "%)", sep = ""), Phylum = paste(sum(nchar(data$tax$Phylum) >
    3), "(", round(sum(nchar(data$tax$Phylum) > 3) / nrow(data$abund) *
    100, digits = 2), "%)", sep = ""), Class = paste(sum(nchar(data$tax$Class) >
    3), "(", round(sum(nchar(data$tax$Class) > 3) / nrow(data$abund) *
    100, digits = 2), "%)", sep = ""), Order = paste(sum(nchar(data$tax$Order) >
    3), "(", round(sum(nchar(data$tax$Order) > 3) / nrow(data$abund) *
    100, digits = 2), "%)", sep = ""), Family = paste(sum(nchar(data$tax$Family) >
    3), "(", round(sum(nchar(data$tax$Family) > 3) / nrow(data$abund) *
    100, digits = 2), "%)", sep = ""), Genus = paste(sum(nchar(data$tax$Genus) >
    3), "(", round(sum(nchar(data$tax$Genus) > 3) / nrow(data$abund) *
    100, digits = 2), "%)", sep = ""), Species = paste(sum(nchar(data$tax$Species) >
    3), "(", round(sum(nchar(data$tax$Species) > 3) / nrow(data$abund) *
    100, digits = 2), "%)", sep = "")), justify = "right")
  cat(
    "\nMetadata variables:", as.character(ncol(data$metadata)),
    "\n", paste(as.character(colnames(data$metadata)), collapse = ", ")
  )
}

#' Rarefaction curve
#'
#' @param data (required) Data list as loaded with amp_load.
#' @param stepsize Step size for the curves. Lower is prettier but takes more time to generate. (default: 1000)
#' @param color_by Color curves by a variable in the metadata.
#'
#' @return A ggplot2 object.
#'
#' @description This function replaces the ampvis2 function amp_rarecurve to fix subsampling labeling bug
#' in vegan
#'
#' @source [utilities.R](../R/utilities.R)
#'
amp_rarecurvefix <- function (data, stepsize = 1000, color_by = NULL) {
  if (class(data) != "ampvis2")
    stop("The provided data is not in ampvis2 format. Use amp_load() to load your data before using ampvis functions. (Or class(data) <- \"ampvis2\", if you know what you are doing.)")
  maxreads <- max(colSums(data$abund))
  if (maxreads < stepsize) {
    stop("\"stepsize\" too high, maximum number of reads in any sample is: ",
         maxreads)
  }
  abund <- t(as.matrix(data[["abund"]]))
  metadata <- data[["metadata"]]
  if (!identical(all.equal(abund, round(abund)), TRUE))
    stop("Function accepts only integers (counts)")
  tot <- rowSums(abund)
  nr <- nrow(abund)
  out <- lapply(seq_len(nr), function(i) {
    if (tot[i] < stepsize) {
      n <- c(1, tot[i])
    } else {
      n <- seq(1, tot[i], by = stepsize)
    }
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    } else {
      n <- c(n[1:(length(n) - 1)], tot[i])
    }
    drop(vegan::rarefy(abund[i, ], n))
  })
  df <- data.frame(Reads = as.numeric(), Species = as.numeric(),
                   SampleID = as.character())
  for (i in 1:length(out)) {
    tsample <- names(attributes(out[[i]])$Subsample[length(out[[i]])])
    tspecies <- unlist(out[[i]])
    treads <- attributes(out[[i]])$Subsample
    tdf <- data.frame(Reads = treads, Species = tspecies,
                      SampleID = tsample)
    df <- rbind.data.frame(df, tdf)
  }
  metadata_col1name <- colnames(metadata)[1]
  colnames(df)[which(colnames(df) == "SampleID")] <- metadata_col1name
  dfm <- merge(metadata, df, by = metadata_col1name)
  p <- ggplot2::ggplot(dfm, ggplot2::aes_string(x = "Reads", y = "Species", group = metadata_col1name,
                              color = color_by)) + ggplot2::geom_line() + ggplot2::theme_classic() +
    ggplot2::xlab("Sequencing depth (reads)") + ggplot2::ylab("Number of observed OTUs")
  return(p)
}


#' biomformat read_biom
#'
#' @param biom_file input biom file name
#'
#' @return biom object
#'
#' @description This function replaces the biomformat function read_biom to deal with reading in
#' crappy hdf5 biom file.
#'
#' @importFrom jsonlite fromJSON
#'
read_biom <- function (biom_file)
{

  errmsg <- paste0("Both attempts to read input file:\n", biom_file,
                   "\n", "either as JSON (BIOM-v1) or HDF5 (BIOM-v2) failed.\n",
                   "Check file path, file name, file itself, then try again.")

  if (!file.exists(biom_file)) {
    stop(paste(biom_file, "does not exist."))
  }
  trash = try(silent = TRUE, expr = {
    x <- fromJSON(biom_file, simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
  })
  if (inherits(trash, "try-error")) {
    tempbiom <- file.path(tempdir(), "temp.biom")
    logoutput(paste("Attempting to convert biom file to", tempbiom))
    trash = try(silent = TRUE, expr = {
      system2("biom", c("convert", "-i", biom_file, "-o", tempbiom, "--to-json", "--header-key", "taxonomy"))
    })

    if (file.exists(tempbiom)) on.exit(file.remove(tempbiom))

    if (inherits(trash, "try-error")) {
      logoutput("file conversion failed")
      stop(errmsg)
    } else {
      trash = try(silent = TRUE, expr = {
        x <- fromJSON(tempbiom, simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
      })

      if (inherits(trash, "try-error")) {
        logoutput("reading from json failed")
        stop(errmsg)
      }
    }
  }

  return(biomformat::biom(x))
}

#' Log base 10 + 1 scale
#'
#' @description Transformation which computes \code{log10(x+1)} scale
#'
#' @name log10scale
#'
#' @return \code{log10p} returns a scales tranformation object
#'
#' @details \code{log10p} is for use with ggplot2 \code{trans} argument in scale function.
#'
log10p_trans <- function() {
  br <- function(x) {
    x <- x+1
    return(round(scales::log_breaks()(x)))
  }
  scales::trans_new("log10p", function(x) log10(x+1) , function(x) (10^(x) - 1), domain = c(0, Inf), format=scales::format_format(digits=4), breaks = br)
}
