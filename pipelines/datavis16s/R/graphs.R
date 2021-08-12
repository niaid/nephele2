#' Subset and rarefy OTU table.
#'
#' @description Subset and/or rarefy OTU table.
#' @source [graphs.R](../R/graphs.R)
#' @param amp  ampvis2 object
#' @param sampdepth  sampling depth.  See details.
#' @param rarefy  rarefy the OTU table in addition to subsetting
#' @param printsummary Logical. print ampvis2 summary of OTU table
#' @param outdir Output directory.  If not null, and samples are removed from amp, the sample names will be output
#' to outdir/samples_being_ignored.txt
#' @param ... other parameters to pass to amp_subset_samples
#'
#' @details \code{sampdepth} will be used to filter out samples with fewer than this number of reads.  If
#' rarefy is TRUE, then it will also be used as the depth at which to subsample using vegan function
#' rrarefy.
#'
#' @return ampvis2 object
#'
#' @importFrom vegan rrarefy
#'
subsetamp <- function(amp, sampdepth = NULL, rarefy=FALSE, printsummary=T, outdir=NULL, ...) {
  ## rarefy
  if (rarefy & !is.null(sampdepth)) {
    cmnd <- "set.seed(500)"
    logoutput(cmnd)
    eval(parse(text=cmnd))
    cmnd <- 'otu <- rrarefy(t(amp$abund), sampdepth)'
    logoutput(cmnd)
    eval(parse(text=cmnd))
    amp$abund <- as.data.frame(t(otu))
  }

  ## subset samples
  if (is.null(sampdepth)) sampdepth = 0
  samples <-  amp$metadata$SampleID
  cmnd <- paste0('amp <- amp_subset_samples(amp, minreads = ', sampdepth,', ...)')
  if (printsummary) logoutput(cmnd)
  eval(parse(text=cmnd))

  ## excluded samples
  if (length(excluded <- setdiff(samples, amp$metadata$SampleID)) > 0) {
    ## output excluded samples to file as well
    if (!is.null(outdir)) {
      excludedpath <- file.path(outdir, "samples_being_ignored.txt")
      logoutput(paste("Saving excluded sample ids to", excludedpath))
      write(excluded, file=excludedpath, ncolumns=1)
    } else {
      writeLines(c("Samples excluded:", excluded))
    }
  }

  ## print ampvis summary
  if (printsummary) {
    print_ampvis2(amp)
    writeLines('')
  }


  return(amp)
}



#' Read in data
#'
#' @param datafile  full path to input data file.  must be either biom file or tab delimited text file.
#' See details.
#' @param mapfile  full path to mapfile.  must contain SampleID, TreatmentGroup, and Description columns
#' @param tsvfile  Logical.  Is datafile a tab-delimited text file? See details.
#' @param mincount  minimum number of reads
#' @return ampvis2 object
#'
#' @details datafile may be either biom file or text file.  If text file, it should have ampvis2 OTU table
#' format \url{https://madsalbertsen.github.io/ampvis2/reference/amp_load.html#the-otu-table}. If the
#' number of reads is less than mincount, the function will give an error, as we cannot make graphs
#' with so few counts.
#'
#' @source [graphs.R](../R/graphs.R)
#'
#' @importFrom utils read.delim
#'
#' @export
#'
readindata <- function(datafile, mapfile, tsvfile=FALSE, mincount=10) {
  #  locale <- Sys.getlocale()
  Sys.setlocale('LC_ALL','C')
  ## mapfile
  logoutput(paste("Reading in map file", mapfile))
  map <- read.delim(mapfile, check.names = FALSE, colClasses = "character",  na.strings = '', comment.char = '')
  colnames(map) <- gsub("^\\#SampleID$", "SampleID", colnames(map))

  if (!all(c("SampleID", "TreatmentGroup", "Description") %in% colnames(map))) {
    stop("Map file missing necessary columns.")
  }

  ## otu table
  logoutput(paste("Reading in OTU file", datafile))
  ## text file
  if ((tsvfile)) {
    cmnd <- "otu <- read.delim(datafile, check.names = FALSE, na.strings = '', row.names = 1)"
    logoutput(cmnd)
    eval(parse(text = cmnd))

    cmnd <- 'tax <- otu[,!names(otu) %in% map$SampleID]'
    logoutput(cmnd)
    eval(parse(text = cmnd))

    cmnd <- 'otu <- otu[, names(otu) %in% map$SampleID, drop=F]'
    logoutput(cmnd)
    eval(parse(text = cmnd))

  } else {

    ## biom file
    cmnd <- 'biom <- read_biom(datafile)'
    logoutput(cmnd)
    eval(parse(text=cmnd))

    cmnd <- 'otu <- as.data.frame(as.matrix(biomformat::biom_data(biom)))'
    logoutput(cmnd)
    eval(parse(text = cmnd))

    cmnd <- 'tax <- biomformat::observation_metadata(biom)'
    logoutput(cmnd)
    eval(parse(text = cmnd))

    ## delete biom structure
    rm(biom)
    gc()

    ## qiime biom file
    if (class(tax) == "list") {
      cmnd <- 'tax <- t(sapply(tax, "[", i=1:7))'
      logoutput(cmnd)
      eval(parse(text = cmnd))
    }

    ## Check for bootstrap columns, i.e. mothur
    if ("bootstrap1" %in% colnames(tax)) {
      bootstrapcols <- grep("bootstrap", colnames(tax))
      tax <- tax[,-bootstrapcols]
    }
  }

  ## Check if tax has 7 columns
  if (ncol(tax) != 7) {
    if (ncol(tax) == 6) {
      tax <- cbind.data.frame(tax, tax[,ncol(tax)])
    } else {
      stop("taxonomy does not have 7 levels.")
    }
  }

  cmnd <- 'otu <- cbind(otu, tax)'
  logoutput(cmnd)
  eval(parse(text = cmnd))


  colnames(otu)[seq.int(to = ncol(otu), length.out = 7)] <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  cmnd <- 'amp <- amp_load(otu, map)'
  logoutput(cmnd)
  eval(parse(text = cmnd))
  print_ampvis2(amp)
  writeLines('')

  if (sum(amp$abund) < mincount) {
    stop(paste0("Not enough counts in OTU table to make any graphs. At least ", mincount, " are needed."))
  }
  #  Sys.setlocale('LC_ALL',locale)
  return(amp)

}


#' Make rarefaction curve graph
#'
#' @param datafile full path to input OTU file (biom or see \link{readindata})
#' @param outdir full path to output directory
#' @param mapfile  full path mapping file
#' @param amp  (Optional) ampvis2 object. may be specified instead of mapfile and datafile
#' @param colors (Optional) color vector - length equal to number of TreatmentGroups in mapfile
#' @param cat Category/column in mapping file by which to color the curves in the graph.
#' (default TreatmentGroup)
#' @param stepsize for rarefaction plotting.
#' @param ... parameters to pass to \code{\link{readindata}}
#'
#' @return Saves rarefaction curve plot to output directory.
#' @export
#'
#' @importFrom plotly ggplotly plotly_data
#' @importFrom ggplot2 scale_color_manual
#'
#' @source [graphs.R](../R/graphs.R)
#'
rarefactioncurve <- function(datafile, outdir, mapfile, amp = NULL, colors=NULL, cat = "TreatmentGroup", stepsize=1000, ...) {

  ## read in data
  if (is.null(amp)) {
    cmnd <- 'amp <- readindata(datafile=datafile, mapfile=mapfile, ...)'
    logoutput(cmnd)
    eval(parse(text = cmnd))
  }
  rarecurve <- amp_rarecurvefix(amp, color_by = cat, stepsize=stepsize) + ggtitle("Rarefaction Curves")

  on.exit(graphics.off())

  ## colors for curves
  if (!is.null(colors)) rarecurve <- rarecurve + scale_color_manual(values = colors)

  ## suppress ggplotly warning to install dev version of ggplot2, as it is out of date
  withCallingHandlers({
    ## plot curves
    rarecurve <- ggplotly(rarecurve, tooltip = c("SampleID", "x", "y"))
   }, message = function(c) {
     if (startsWith(conditionMessage(c), "We recommend that you use the dev version of ggplot2"))
       invokeRestart("muffleMessage")
   })

  ## make table of data for plotly export
  df <- plotly_data(rarecurve)
  df_new <- split(df, df$SampleID)
  tg <- which(colnames(df) == "TreatmentGroup")
  desc <-  which(colnames(df) == "Description")
  n <- ncol(df)

  df_new <- lapply(df_new, function(x) {x <- data.frame(x); y <- unlist(c(rep(NA, tg - 1), x[1,tg:desc], rep(NA,(n - desc)))); names(y) <- colnames(x); rbind(x,y)   })
  df_new <- do.call(rbind.data.frame, df_new)


  ## save to html and txt file
  plotlyGrid(rarecurve, file.path(outdir, "rarecurve.html"), data = df_new)
  logoutput(paste0('Saving rarefaction curve table to ', file.path(outdir, 'rarecurve.txt') ))
  write.table(df_new, file.path(outdir, 'rarecurve.txt'), quote = FALSE, sep = '\t', row.names = FALSE, na = "")

  return(as.integer(0))

}


#' PCoA plots
#'
#' @param datafile full path to input OTU file (biom or see \link{readindata})
#' @param outdir  full path to output directory
#' @param mapfile  full path to map file
#' @param amp  ampvis2 object. may be specified instead of mapfile and datafile
#' @param sampdepth  sampling depth
#' @param distm  distance measure for PCoA.  any that are supported by
#' \href{https://madsalbertsen.github.io/ampvis2/reference/amp_ordinate.html}{amp_ordinate} except for unifrac, wunifrac, and none.
#' @param filter_species Remove low abundant OTU's across all samples below this threshold in percent.
#' Setting this to 0 may drastically increase computation time.
#' @param rarefy Logical. Rarefy the OTU table if sampdepth is specified.
#' @param colors (Optional) color vector - length equal to number of TreatmentGroups in mapfile
#' @param filesuffix (Optional) suffix for output filename
#' @param ... parameters to pass to  \code{\link{readindata}}
#'
#' @return Saves pcoa plots to outdir.
#'
#' @importFrom ggplot2 scale_color_manual scale_fill_manual ggtitle
#'
#' @source [graphs.R](../R/graphs.R)
#'
#' @export
#'
pcoaplot <- function(datafile, outdir, mapfile, amp=NULL, sampdepth = NULL, distm="binomial", filter_species=0.1, rarefy=FALSE, colors=NULL, filesuffix=NULL, ...) {

  ## read in data
  if (is.null(amp)) {
    cmnd <- 'amp <- readindata(datafile=datafile, mapfile=mapfile, ...)'
    logoutput(cmnd)
    eval(parse(text = cmnd))
  }

  ## remove samples below sampdepth and rarefy, if necessary
  if (!is.null(sampdepth)) {
    cmnd <- paste0('amp <- subsetamp(amp, sampdepth = ', sampdepth,', rarefy=',rarefy, ')')
    logoutput(cmnd)
    eval(parse(text=cmnd))
  }

  on.exit(graphics.off())

  ## plot PCoA
  cmnd <- paste0('pcoa <- amp_ordinate(amp, filter_species =', filter_species, ',type="PCOA", distmeasure ="', distm, '",sample_color_by = "TreatmentGroup", detailed_output = TRUE, transform="none")')

  logoutput(cmnd)
  eval(parse(text = cmnd))
  if (!is.null(colors)) pcoa$plot <- pcoa$plot + scale_color_manual(values = colors) + scale_fill_manual(values=colors) + ggtitle(paste("PCoA with", distm, "distance"))

  ## save to file
  outfile <- file.path(outdir, paste0("pcoa_", distm, filesuffix, ".html"))
  plotlyGrid(pcoa$plot, outfile, data = pcoa$dsites)
  tabletsv <- gsub('.html$', '.txt', outfile)
  logoutput(paste0('Saving ', distm, ' PCoA table to ', tabletsv))
  write.table(pcoa$dsites, tabletsv, quote=FALSE, sep='\t', row.names=FALSE, na="")

  return(as.integer(0))

}


#' Morpheus heatmap
#'
#' @description Creates heatmaps using Morpheus R API \url{https://software.broadinstitute.org/morpheus/}.  The heatmaps are made
#' using relative abundances.
#'
#' @param datafile  full path to input OTU file (biom or see \link{readindata})
#' @param outdir  full path to output directory
#' @param mapfile full path to mapping file
#' @param amp  (Optional) ampvis2 object. may be specified instead of mapfile and datafile
#' @param sampdepth sampling depth
#' @param rarefy Logical. Rarefy the OTU table if sampdepth is specified.
#' @param filter_level minimum abundance to show in the heatmap
#' @param taxlevel vector of taxonomic levels to graph.  must be subset of
#' c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "seq").  See Details.
#' @param colors  (Optional) color vector - length equal to number of TreatmentGroups in mapfile
#' @param rowAnnotations (Optional) Row annotations to be used in addition to taxonomy.
#' @param force Force "seq" level heatmap to be made even if number of seqs is greater than 2000. '
#' See Details.
#' @param filesuffix (Optional) suffix for output filename
#' @param ...  parameters to pass to  \code{\link{readindata}}
#'
#' @return  Saves heatmaps to outdir.
#' @export
#'
#' @details For the \code{taxlevel} parameter, each level is made into a separate heatmap.  "seq" makes
#'  the heatmap with no collapsing of taxonomic levels if there are fewer than 2000 ASVs/OTUs.  Otherwise,
#'  Species level is made instead.
#'
#' @importFrom morpheus morpheus
#' @importFrom grDevices colorRampPalette
#' @importFrom htmltools tags
#'
#' @source [graphs.R](../R/graphs.R)
#'
#' @examples
#' \dontrun{
#'  morphheatmap(datafile="OTU_table.txt", outdir="outputs/graphs", mapfile="mapfile.txt",
#'  sampdepth = 25000, taxlevel = c("Family", "seq"), tsvfile=TRUE)
#' }
#'
#'
morphheatmap <- function(datafile, outdir, mapfile, amp = NULL, sampdepth = NULL, rarefy=FALSE, filter_level = NULL, taxlevel=c("seq"), colors = NULL, rowAnnotations=NULL, force=FALSE, filesuffix=NULL, ...) {

  ## read in data
  if (is.null(amp)) {
    cmnd <- 'amp <- readindata(datafile=datafile, mapfile=mapfile, ...)'
    logoutput(cmnd)
    eval(parse(text = cmnd))
  }


  ## filter low abundant taxa, if desired
  if (!is.null(filter_level)) {
    logoutput(paste('Filter taxa below', filter_level, 'counts/abundance.'))
    if (filter_level < 1) {
      cmnd <- paste0('amp <- filterlowabund(amp, level = ', filter_level,')')
      logoutput(cmnd)
      eval(parse(text = cmnd))
    } else {
      cmnd <- paste0('amp <- filterlowabund(amp, level = ', filter_level,', abs=T)')
      logoutput(cmnd)
      eval(parse(text = cmnd))
    }
  }


  ## normalize data
  logoutput("Calculate relative abundance.")
  cmnd <- paste0("amp <- subsetamp(amp, sampdepth = ", deparse(sampdepth) ,", rarefy=", rarefy, ", normalise = TRUE, printsummary = FALSE)")
  logoutput(cmnd)
  eval(parse(text=cmnd))


  if (!all(nrow(amp$abund) > 1, ncol(amp$abund) > 1)) {
    stop("OTU table must be at least 2x2 for heatmap.")
  }

  ## Annotations
  tg <- which(colnames(amp$metadata) == "TreatmentGroup")
  desc <-  which(colnames(amp$metadata) == "Description") - 1
  columns <- lapply(tg:desc, function(x) { list(field=colnames(amp$metadata)[x], highlightMatchingValues=TRUE, display=list('color'))  } )
  columns <- append(list(list(field='id', display=list('text'))), columns)
  rows <- list(list(field='id', display=list('text')))

  ## Heatmap colors
  hmapcolors <- rev(colorRampPalette( RColorBrewer::brewer.pal(11, "RdYlBu"), bias=0.5)(100))

  ## heatmap function
  makeheatmap <- function(tl, amp) {

    ## create matrix at higher taxonomic level
    if (tl != "seq") {
      amptax <- highertax(amp, taxlevel=tl)
    } else {
      amptax <- amp
      amptax$tax <- shortnames(amptax$tax)
    }

    ## row and column names for matrix
    if (tl == "seq") {
      sn <- paste(amptax$tax$OTU, amptax$tax$Species)
    } else {
      sn <- amptax$tax[,ncol(amptax$tax)]
    }
    mat <- amptax$abund
    row.names(mat) <- sn
    mat <- mat[,amptax$metadata$SampleID]
    if (!is.null(rowAnnotations)) {
      amptax$tax <- cbind.data.frame(amptax$tax, rowAnnotations[match(row.names(amptax$tax), rowAnnotations$taxa),])
    }

    ## log scale for colors
    mm <- max(amptax$abund)
    minm <- min(amptax$abund[which(amptax$abund > 0, arr.ind = T)])
    values <-  c(0,expm1(seq(log1p(minm), log1p(100), length.out = 99)))
    w <- which(values > 10)
    values[w] <- round(values[w])


    ## make morpheus heatmap
    cmnd <- 'heatmap <- morpheus(mat, columns=columns, columnAnnotations = amptax$metadata, columnColorModel = list(type=as.list(colors)), colorScheme = list(scalingMode="fixed", values=values, colors=hmapcolors, stepped=FALSE), rowAnnotations = amptax$tax, rows = rows, dendrogram="none")'
    logoutput(cmnd)
    eval(parse(text = cmnd))

    ## Save html file
    outdir <- tools::file_path_as_absolute(outdir)
    outfile <- file.path(outdir, paste0(tl, "_heatmap",filesuffix, ".html"))
    logoutput(paste("Saving plot to", outfile))
    heatmap$width = '100%'
    heatmap$height = '90%'
    tt <- tags$div(heatmap, style = "position: absolute; top: 10px; right: 40px; bottom: 40px; left: 40px;")
    save_fillhtml(tt, file = outfile, bodystyle = 'height:100%; width:100%;overflow:hidden;')

  }

  ## make heatmap at different taxonomic levels
  for (t in taxlevel) {
    if (nrow(amp$abund) > 2000 & t == "seq" & !force) {
      logoutput("Number of sequence variants > 2000. Making heatmap at species/lowest assigned taxonomic level.")
      t = "Species"
    }

    cmnd <- paste0('makeheatmap("', t, '", amp)')
    logoutput(cmnd)
    if (inherits(try(eval(parse(text = cmnd))), "try-error")) {
      stop(paste("Heatmap at ", t, " level failed."))
    }
  }

  return(as.integer(0))

}


#' Alpha diversity boxplot
#'
#' Plots plotly boxplot of shannon diversity and Chao species richness.  If sampling depth is NULL,
#' rarefies OTU table to the minimum readcount of any sample.  If this is low, then the plot will fail.
#'
#' @param datafile full path to input OTU file
#' @param outdir  full path to output directory
#' @param mapfile  full path to map file
#' @param amp  ampvis2 object. may be specified instead of mapfile and datafile
#' @param sampdepth  sampling depth.  see details.
#' @param colors colors to use for plots
#' @param cats categories/columns in mapping file to use as groups.  If NULL (default), will use
#' all columns starting with TreatmentGroup to (but not including) Description
#' @param filesuffix (Optional) suffix for output filename
#' @param ... other parameters to pass to \link{readindata}
#'
#' @return Save alpha diversity boxplots to outdir.
#' @export
#'
#' @details If \code{sampdepth} is NULL, the sampling depth is set to the size of the smallest
#' sample.
#'
#' @source [graphs.R](../R/graphs.R)
#'
#' @importFrom plotly subplot ggplotly
#' @importFrom ggplot2 scale_color_manual element_text
#'
adivboxplot <- function(datafile, outdir, mapfile, amp=NULL, sampdepth = NULL, colors = NULL, cats = NULL, filesuffix=NULL, ...) {

  ## read in data
  if (is.null(amp)) {
    cmnd <- 'amp <- readindata(datafile=datafile, mapfile=mapfile, ...)'
    logoutput(cmnd)
    eval(parse(text = cmnd))
  }

  ## rarefy to minimum sample count if sampling depth is not specified
  if (is.null(sampdepth)) {
    logoutput('Calculate number of counts to rarefy table.')
    cmnd <- 'sampdepth <-  min(colSums(amp$abund))'
    logoutput(cmnd)
    eval(parse(text=cmnd))
  }

  ## compute alpha diversity and species richness
  cmnd <- paste0('alphadiv <- amp_alphadiv(amp, measure="shannon", richness = TRUE, rarefy = ', sampdepth, ')')
  logoutput(cmnd)
  eval(parse(text = cmnd))

  ## save tabular output
  logoutput(paste0('Saving alpha diversity table to ', file.path(outdir, 'alphadiv.txt') ))
  write.table(alphadiv, file.path(outdir, 'alphadiv.txt'), quote = FALSE, sep = '\t', row.names = FALSE, na = "")

  divplots <- function(adiv, col, colors) {
    if (col == "TreatmentGroup") {
      lc <- unname(colors[1:length(unique(adiv$TreatmentGroup))])
    } else {
      lc <- NULL
    }

    gshan <- ggplot2::ggplot(adiv, ggplot2::aes_string(x=col, y="Shannon", color = col, text = "SampleID")) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_jitter(width = 0.2, height = 0, size=1) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "none", axis.text.x = element_text(angle=60)) + ggtitle("species alpha diversity")
    if (!is.null(lc)) gshan <- gshan + scale_color_manual(values=lc)
    shannon <- ggplotly(gshan, tooltip = c("text", "x", "y"))

    gchao1 <- ggplot2::ggplot(adiv, ggplot2::aes_string(x=col, y="Chao1", color = col, text = "SampleID")) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_jitter(width = 0.2, height = 0, size=1) + ggplot2::theme_bw() + ggplot2::theme(legend.position = "none", axis.text.x = element_text(angle=60))
    if (!is.null(lc)) gchao1 <- gchao1 + scale_color_manual(values=lc)
    chao1 <- ggplotly(gchao1, tooltip = c("text", "x", "y"))

    return(list(shannon, chao1))
  }

  if (is.null(cats)) {
    ## make plots for each category (between treatmentgroup and description)
    tg <- which(colnames(amp$metadata) == "TreatmentGroup")
    desc <- which(colnames(amp$metadata) == "Description") - 1
    cats <- colnames(amp$metadata)[tg:desc]
  }

  ## style widgets 2 across
  plotlyplots <- unlist(lapply(cats, function(x) divplots(alphadiv, x, colors)), recursive = FALSE)

  ### manually remove outliers https://community.plotly.com/t/ggplotly-ignoring-geom-boxplot-outlier-parameters/2247
  plotlyplots <- lapply(plotlyplots, function(p) {
    p$x$data[1:length(p$x$data)/2] <- lapply(p$x$data[1:length(p$x$data)/2],FUN=function(y){y$marker=list(opacity=0); return(y)})
    return(p)
  })

  divwidget <- plotly::subplot(plotlyplots, nrows = length(plotlyplots)/2, titleX=TRUE, titleY = TRUE, margin = 0.06)
  if (length(cats) > 2) divwidget$height <- paste0(400*(length(cats)), "px")

  ## save html file
  plotlyGrid(divwidget, file.path(outdir, "alphadiv.html"), data = alphadiv)
  return(as.integer(0))

}

#' Pipeline function
#'
#' @description Make all 4 types of graphs

#' @param datafile full path to input OTU file (biom or txt file see \link{readindata} for format)
#' @param outdir  full path to output directory
#' @param mapfile  full path to map file
#' @param sampdepth  sampling depth.  default: 10000
#' @param ... other parameters to pass to \link{readindata}
#'
#' @return graphs are saved to outdir.  See [user doc](../doc/user_doc.md).
#'
#' This value is used to remove samples before for alpha diversity and PCoA plots.
#' Also, to rarefy OTU table for the alpha diversity and Bray-Curtis distance PCoA.
#'
#' @importFrom stats median
#'
#' @source [graphs.R](../R/graphs.R)
#'
#' @export
#'
allgraphs <- function(datafile, outdir, mapfile, sampdepth = 10000, ...) {

  ## Set initial return value to 0.
  retvalue <- as.integer(0)

  ## Read in abundance data and mapfile
  amp <- readindata(datafile=datafile, mapfile=mapfile, ...)

  ## Choose colors
  amp$metadata <-  amp$metadata[order(amp$metadata$TreatmentGroup),]
  numtg <- length(unique(amp$metadata$TreatmentGroup))
  allcols <- c(RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"))[c(1:6,9,11:13,16, 15,14,18:20)]
  coln <- length(allcols)
  if (numtg <= coln){
    allcols <- allcols[1:numtg]
    names(allcols) <- unique(amp$metadata$TreatmentGroup)
  } else {
    allcols <- NULL
  }

  ## Rarefaction curve
  logoutput('Rarefaction curve', 1)
  cmnd <- 'rarefactioncurve(outdir = outdir, amp = amp, colors = allcols)'
  logoutput(cmnd)
  if (inherits(try(eval(parse(text=cmnd))), "try-error")) retvalue <- as.integer(1)

  ## Heatmap
  logoutput('Relative abundance heatmaps', 1)
  if (!all(nrow(amp$abund) > 1, ncol(amp$abund) > 1)) {
    logoutput(paste("OTU table needs to be at least 2x2 for heatmaps.  Only", nrow(amp$metadata), "samples and", nrow(amp$abund) ,"sequence variants are in the dataset.  Heatmaps and diversity plots will not be made."), 1, type='WARNING')
    return(retvalue)
  }


  cmnd <- 'morphheatmap(outdir = outdir, amp = amp, colors=allcols, filter_level = 5)'
  logoutput(cmnd)
  if (inherits(try(eval(parse(text=cmnd))), "try-error")) retvalue <- as.integer(1)

  ## Filter out low count samples
  logoutput(paste('Sampling depth:', sampdepth), 1)
  cs <- colSums(amp$abund)
  if (max(cs) < sampdepth) {
    logoutput(paste("The counts for all samples are below the sampling depth of", sampdepth, ", so diversity and PCoA plots will not be made."), 1, type="WARNING")
    return(retvalue)
  }

  logoutput(paste0('Filter samples below ', sampdepth, ' counts.'))
  ampsub <- subsetamp(amp, sampdepth=sampdepth, outdir = outdir)

  if (nrow(ampsub$metadata) < 3) {
    logoutput(paste("Alpha diversity and PCoA plots will not be made, as they require at least 3 samples.  Only", nrow(ampsub$metadata), "remain after filtering."), 1, type='WARNING')
    return(retvalue)
  }


  ## binomial PCoA
  logoutput('PCoA plot with binomial distance', 1)
  cmnd <- paste0('pcoaplot(outdir = outdir, amp = ampsub, distm = "binomial", colors = allcols)')
  logoutput(cmnd)
  if (inherits(try(eval(parse(text=cmnd))), "try-error")) retvalue <- as.integer(1)



  ## Rarefy table
  logoutput(paste('Rarefying OTU Table to ', sampdepth, 'reads.'), 1)
  amprare <- subsetamp(ampsub, sampdepth = sampdepth, rarefy = TRUE, normalise=FALSE, printsummary = T)
  rareotutablefile <-  file.path(outdir, paste0('rarefied_OTU_table_', sampdepth, '.txt'))
  logoutput(paste('Saving rarefied OTU Table to ',   rareotutablefile ))
  rareotutable <- cbind.data.frame(amprare$abund, amprare$tax)
  rareotutable$OTU <- NULL
  write.table(rareotutable, rareotutablefile , quote = FALSE, sep = '\t', col.names = NA, na = "")


  logoutput("Making heatmap from rarefied counts.")
  cmnd <- 'morphheatmap(outdir = outdir, amp = amprare, colors=allcols, filter_level = 5, filesuffix = "_rarefied")'
  logoutput(cmnd)
  if (inherits(try(eval(parse(text=cmnd))), "try-error")) retvalue <- as.integer(1)



  ## bray-curtis PCoA
  logoutput(paste('Normalizing rarefied OTU table to 100 for Bray-Curtis distance.'), 1)
  ampbc <- subsetamp(amprare, sampdepth = sampdepth, normalise=TRUE, printsummary = F)

  cmnd <- 'pcoaplot(outdir = outdir, amp = ampbc, distm = "bray", colors = allcols, filesuffix="_rarefied")'
  logoutput(cmnd)
  if (inherits(try(eval(parse(text=cmnd))), "try-error")) retvalue <- as.integer(1)


  ## Alpha diversity
  logoutput('Alpha diversity boxplot', 1)
  cmnd <- 'adivboxplot(outdir = outdir, amp = amprare, sampdepth = sampdepth, colors = allcols)'
  logoutput(cmnd)
  if (inherits(try(eval(parse(text=cmnd))), "try-error")) retvalue <- as.integer(1)


  return(retvalue)

}

#' Wrapper for any graph function
#'
#' @description This is a wrapper for any of the graph functions meant to be called using rpy2 in python.
#'
#' @param datafile full path to input OTU file (biom or txt, see \link{readindata} for format of txt file)
#' @param outdir  output directory for graphs
#' @param mapfile full path to map file
#' @param FUN character string. name of function you would like to run. can be actual
#' function object if run from R
#' @param logfilename logfilename
#' @param info print sessionInfo to logfile
#' @param tsvfile Is datafile a tab-delimited text file? Default FALSE
#' @param ...  parameters needed to pass to FUN
#'
#' @return Returns 0 if FUN succeeds and stops on error.  In rpy2, it will throw
#' rpy2.rinterface.RRuntimeError.
#'
#' @export
#'
#' @source [graphs.R](../R/graphs.R)
#'
#' @importFrom utils capture.output sessionInfo
#'
#' @examples
#'
#' \dontrun{
#'
#' # example with no optional arguments for running allgraphs
#' trygraphwrapper("/path/to/outputs/out.biom", "/path/to/outputs/",
#' "/path/to/inputs/mapfile.txt", 'allgraphs')
#'
#' # example with sampdepth argument for running allgraphs
#' trygraphwrapper("/path/to/outputs/out.biom", "/path/to/outputs/",
#' "/path/to/inputs/mapfile.txt", 'allgraphs', sampdepth=30000)
#'
#'
#' # example with optional argument sampdepth and tsv file
#' trygraphwrapper("/path/to/outputs/OTU_table.txt", "/path/to/outputs/",
#' "/path/to/inputs/mapfile.txt", 'allgraphs', sampdepth = 30000, tsvfile=TRUE)
#'
#' # example of making heatmap with optional arguments
#' trygraphwrapper("/path/to/outputs/taxa_species.biom", "/path/to/outputs",
#' "/path/to/inputs/mapfile.txt", 'morphheatmap', sampdepth = 30000, filter_level=0.01,
#' taxlevel=c("Family", "seq"))
#' }
#'
#'
trygraphwrapper <- function(datafile, outdir, mapfile, FUN, logfilename="logfile.txt", info = TRUE, tsvfile=FALSE, ... ) {

  ## set error handling options here, since rpy2 does not allow setting globally
  ## see http://ai-bcbbsptprd01.niaid.nih.gov:8080/browse/NPHL-769
  options(stringsAsFactors = FALSE, scipen = 999, warn=1, show.error.locations= TRUE, error = function() traceback(2), digits.secs = 3)

  ## open log file
  logfile <- file(logfilename, open = "at")
  sink(file = logfile, type="output")
  sink(file = logfile, type= "message")

  on.exit(closeAllConnections())
  ## create output directory
  outdir <- file.path(outdir, "graphs")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


  ## print sessionInfo
  if (info) writeLines(capture.output(sessionInfo()))

  ## make function command
  functionstring <- deparse(substitute(FUN))
  cmnd <- paste0(functionstring, '(datafile="', datafile, '", outdir="', outdir, '", mapfile="', mapfile, '",', 'tsvfile=',tsvfile, ', ...)')
  logoutput(cmnd, 1)
  FUN <- match.fun(FUN)


  ## run command
  tryCatch(retvalue <- eval(parse(text=cmnd)),
           finally = if (!exists("retvalue")) {
             logoutput(paste(functionstring, "completed with some errors."), bline=1, aline=1, type='ERROR')
           } else if (retvalue == 1) {
             logoutput(paste(functionstring, "completed with some warnings."), bline=1, aline=1, type='WARNING')
           } else {
             logoutput(paste(functionstring,"complete."), bline=1, aline=1)
           })

  return(as.integer(0))

}
