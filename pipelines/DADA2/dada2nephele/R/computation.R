#' write log output
#'
#' @description Prints time along with log message.
#'
#' @param c String. Log message/command to print.
#' @param bline Number of blank lines to precede output.
#' @param aline Number of blank lines to follow output.
#' @param type  String.  Must be one of "WARNING", or "ERROR" or NULL.
#'
#' @source [computation.R](../R/computation.R#L10)
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

#' wrap command in tryCatch
#'
#' @description (optionally) log cmd to output and evaluate cmd in the parent environment, catching
#' and parsing errors.
#'
#' @param cmd Character. Command string.
#' @param step (Optional). Step name to pass to error.  Default will use cmd.
#' @param bline Number of blank lines to precede log output; parameter for \code{\link{logoutput}}
#' @param aline Number of blank lines to follow log output; parameter for \code{\link{logoutput}}
#' @param w2e warning message to escalate to error.
#' @param log log command to file
#'
#' @return cmd will be evaluated in the parent environment, so return values in cmd will be there.
#'
#' @source [computation.R](../R/computation.R#L40)
#'
run_cmd <- function(cmd, step=NULL, bline=0, aline=0, w2e=NA, log=T) {
   error_function <- function(e, logstack=T) {
       if (logstack)
        try({
          sc <- sys.calls()
          realstack <- unique(as.character(c(sc[8:(length(sc)-2)], e$call)))
          logoutput(paste0("Sys.calls:\n", paste(length(realstack):1, realstack, sep=": ", collapse='\n')), type="ERROR")
          }, silent=TRUE)
    if (!is.null(step)) {
      callstring <- paste(deparse(e$call), collapse= " ", sep=" ")
      callstring <- gsub("\\s+", " ", callstring)
      fakestack <- paste0("\nCall: ", callstring, ", Pipeline Step: ", step, ', Pipeline: dada2compute\n\n')
      e$message <- paste0('\n ', e$message)
      msglength <- getOption("warning.length") - nchar(fakestack)
      if (nchar(e$message) > msglength) e$message <- substr(e$message, 1, msglength)
      e$message <- paste0(e$message, fakestack)
      e$call <- NULL
    }
    stop(e)
  }
  if (log) logoutput(cmd, bline=bline, aline=aline)
  if(is.null(step)) step <- cmd
  if (is.na(w2e))
      withCallingHandlers(eval.parent(parse(text=cmd)), error = error_function)
  else
      withCallingHandlers(eval.parent(parse(text=cmd)), error = error_function, warning = function(w) {
          if (grepl(w2e, w$message))
              error_function(w, logstack=F)
      })
}


#' Check files after filterAndTrim
#'
#' @param A  mapping data.frame
#' @param filt.dir filtered data directory
#' @param trimlist list of vectors, R1 and R2 of trimmed files
#'
#' @return list of `A`, `trimr1`, and maybe `trimr2` with missing files/rows removed.
#' Stops on error if no trimmed files exist.
#'
#' @source [computation.R](../R/computation.R#L83)
#' @md
#'
checktrimfiles <- function(A, filt.dir, trimlist) {

  nameslist <- names(trimlist)
  trimlist <- lapply(trimlist, function(x) file.path(filt.dir, x))
  missing.index <- lapply(trimlist, function (x) {
    mf <- which(!file.exists(x))
    return(mf)
  })

  misfiles <- unlist(lapply(nameslist, function(x) trimlist[[x]][missing.index[[x]]] ))
  missing.index <- unique(unlist(missing.index))
  if (length(missing.index) == nrow(A))
    stop("No trimmed files were found.  (For paired end, at least one trimmed FASTQ is missing for each sample)")

  if (length(misfiles) > 0) {
    ## warnfe <- paste(c("The following trimmed files were not found:", misfiles), collapse="\n")
    ## warning(warnfe)
    A <- A[-missing.index,]
    trimlist <- lapply(trimlist, function(x) x[-missing.index])
  }

  return(list(A=A, trimlist = trimlist))
}

#' Combine function foreach loop
#'
#' @description For a \code{\link[foreach]{foreach}} loop that returns a list of n items for
#' each iteration, this function combines all iterations into n different lists - one for each
#' item returned.
#'
#' @param x list of n items to append onto
#' @param ... individual lists of n items to append onto n lists in x
#'
#' @details Use this as the value of the `.combine` parameter in foreach.  You will
#' need to also initiate the list of n lists for the `.init` foreach parameter.
#'
#' @return list of lists named according to names of list from single iteration.
#'
#' @examples
#'
#' \dontrun{
#'   `%do%` <- foreach::`%do%`
#'   oper <- foreach::foreach(i=1:5, .combine='comb', .multicombine=TRUE,
#'   .init=list(list(), list())) %do% {
#'     list(i*2, i*3)
#'   }
#'
#'   oper[[1]]
#' }
#'
#'
#'
comb <- function(x, ...) {

  cc <- lapply(seq_along(x),
               function(i) {c(x[[i]], lapply(list(...), function(y) y[[i]])) })
  names(cc) <- names(x)
  return(cc)
}

#' Get number of reads
#'
#' @description get number of reads for each sample from a dada2 object (such as output of
#' \code{\link[dada2]{mergePairs}}).
#'
#' @param x dada2 object
#'
#' @return Integer vector
#'
getN <- function(x) {
  sum(getUniques(x))
}

#' Assign taxonomy using DECIPHER
#'
#' @param refdb reference database .RData file
#' @param seqtab sequence table made by DADA2
#' @param nthread number of processors to use
#'
#' @return taxa table in DADA2 taxa format
#'
#' @importFrom Biostrings DNAStringSet
#' @importFrom DECIPHER IdTaxa
#'
#'
decipher_assign <- function(refdb, seqtab, nthread) {
  if (grepl("homd", refdb, ignore.case=T)) {
    stop('DECIPHER cannot use HOMD for taxonomic assignment.')
  } else if (!grepl(".RData", refdb)) {
    logoutput(paste('DECIPHER may have problems with', refdb, 'as it must load database from .RData formatted file.'), type="WARNING")
  }
  run_cmd('dna <- DNAStringSet(getSequences(seqtab))', step= "DECIPHER BioStrings::DNAStringSet", bline=1)
  logoutput(paste('Loading refdb', refdb))
  load(refdb)

  ## decipher uses different parallel syntax than dada2
  if (is.numeric(nthread)) {
    nthread=nthread
  } else if (isTRUE(nthread)) {
    nthread=NULL
  } else if (isFALSE(nthread)) {
    nthread=1
  }

  run_cmd('ids <- IdTaxa(dna, trainingSet, strand="both", processors=nthread, verbose=F)', 'DECIPHER::IdTaxa')

  logoutput('reformat Taxa class "ids" to dada2 matrix format "taxid"')
  run_cmd('ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")')
  run_cmd('taxid <- t(sapply(ids, function(x) { m <- match(ranks, x$rank); taxa <- x$taxon[m];  taxa[startsWith(taxa, "unclassified_")] <- NA; taxa }))')
  run_cmd('colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab)', 'DECIPHER::IdTaxa')
  run_cmd('taxa <- taxid')
  return(taxa)
}

#' run dada2 pipeline
#'
#' @rdname dada2compute
#' @md
#'
#' @description \code{dada2compute} will run the whole pipeline from within R.
#' See  \code{\link{.onLoad}} examples for argument default values
#' or \code{getOption('dparams')}.
#'
#' @param datadir Directory containing FASTQ files.
#' @param outdir Output directory.
#' @param mapfile Mapping file.
#' @param refdb Path to the reference database.
#' @param refdb_species Path to the species reference database.
#' @param nthread Number of processors for parallel steps.
#' @param data_type 'SE' or 'PE'
#' @param trimLeft (Optional).  In [dada2::filterAndTrim], the number of nucleotides to remove
#' from the start of each read, forward and reverse.  The values should be chosen based on
#' the lengths of primers used for sequencing. (Whole number vector of length 2, for forward
#' and reverse).
#' @param trimOverhang (Optional). After merging paired end reads, trim sequence which overhangs
#' the start of each read. If amplicons are shorter than read length, e.g. 16S V4 region, we
#' suggest setting this to True.  (Logical. Default: False).
#' @param chimera Should chimeric sequences be removed?  If primers are not
#' @param maxEE Integer. Remove reads with less than maxEE expected errors in [dada2::filterAndTrim]
#' @param truncQ Integer. Truncate reads at first qual below truncQ in [dada2::filterAndTrim]
#' @param minBoot Integer.  Bootstrap value for [dada2::assignTaxonomy]
#' @param truncLen Integer.  Truncation length for reads in [dada2::filterAndTrim].  Can be vector of length 1 or 2
#' @param maxMismatch Integer. Maximum number of allowed mismatches in [dada2::mergePairs]
#' @param justConcatenate Logical. Should PE reads be concatenated instead of overlapped and merged.
#' @param taxmethod method for taxonomic assignment.  either idtaxa or rdp.
#' @param band_size Integer.  Band size for [dada2::dada], [dada2::setDadaOpt]
#' @param homopolymer_gap_penalty Integer. for [dada2::dada] & [dada2::setDadaOpt]. cost of gaps in homopolymer
#' regions (>=3 repeated bases)
#' @param plotquality Logical.  Should quality plots be made?
#'
#' @importFrom foreach foreach "%do%"
#' @importFrom ggplot2 ggsave
#' @importFrom ShortRead writeFasta
#' @importFrom utils capture.output read.delim sessionInfo write.table lsf.str
#' @importFrom stringr str_match
#' @importFrom methods show
#'
#' @source [computation.R](../R/computation.R#L243)
#'
#' @export
#'
dada2compute <- function(datadir,
                         outdir,
                         mapfile,
                         refdb = getOption('dparams')$refdb,
                         refdb_species = getOption('dparams')$refdb_species,
                         nthread = TRUE,
                         chimera = FALSE,
                         trimLeft = getOption("dparams")$trimLeft,
                         trimOverhang = getOption("dparams")$trimOverhang,
                         data_type = "PE",
                         maxEE = getOption('dparams')$maxEE,
                         truncQ = getOption('dparams')$truncQ,
                         minBoot = getOption('dparams')$minBoot,
                         truncLen = getOption('dparams')$truncLen,
                         maxMismatch = getOption('dparams')$maxMismatch,
                         justConcatenate=getOption('dparams')$justConcatenate,
                         taxmethod = getOption('dparams')$taxmethod,
                         band_size = NULL,
                         homopolymer_gap_penalty = NULL,
                         plotquality=T) {


  ## Output session Info
  writeLines(capture.output(sessionInfo()))

  ## Output database information
  writeLines(c("","Taxonomic Reference Database", refdb, refdb_species))

  ## Read in mapping file
  logoutput(paste("Reading in map file ", mapfile), 1)
  A <- read.delim(mapfile, check.names = FALSE, colClasses = "character", na.strings = '', comment.char = '')
  colnames(A) <- gsub("^\\#SampleID$", "SampleID", colnames(A))
  A <- subset(A, !is.na(SampleID))

  ## Parameters
  dparams <- getOption("dparams")


  # dada options ------------------------------------------------------------
  ## Set dada parameters for Ion Torrent data
  if (!is.null(band_size)) setDadaOpt(BAND_SIZE = as.numeric(band_size))
  if (!is.null(homopolymer_gap_penalty)) setDadaOpt(HOMOPOLYMER_GAP_PENALTY = as.numeric(homopolymer_gap_penalty))
  ## print default dada2 options
  logoutput("Printing dada algorithm options.")
  print(format(getDadaOpt(), scientific = 0), quote = F)


  ## Output directory for filtered and trimmed data
  outdir <- normalizePath(outdir)
  filt.dir <- file.path(outdir, "filtered_data")
  ## Output directory for intermediate files
  interm.dir <- file.path(outdir, "intermediate_files")
  dir.create(interm.dir, showWarnings = FALSE)

  if (data_type == "PE") {
    logoutput('Paired End', 1)
    readslist <- list(R1=A$ForwardFastqFile, R2=A$ReverseFastqFile)
  } else {
    logoutput('Single End', 1)
    readslist <- list(R1=A$ForwardFastqFile)
  }

  nameslist <- names(readslist)


# Plot quality profile ----------------------------------------------------
  if (plotquality) {
      if (length(readslist$R1) >= 100) {
          logoutput("Number of samples is >= 100, so we will plot the quality profiles in aggregate", 1)
          cmd <- "pqp <- lapply(readslist, FUN = function(x) { ppp <- plotQualityProfile(file.path(datadir, x), aggregate=T); ppp$facet$params$ncol <- 4; ppp })"
          qplotheight <- 8
      } else {
          writeLines("")
          cmd <- "pqp <- lapply(readslist, FUN = function(x) { ppp <- plotQualityProfile(file.path(datadir, x)); ppp$facet$params$ncol <- 4; ppp })"
          qplotheight <- min(2.5*ceiling(length(readslist[[1]])/4), 49)
      }
      logoutput(cmd)


      ## Wrap quality profile plotting in try statement, so entire pipeline does not fail, if this step fails
      checkpqploterror <- try( {
          eval(parse(text=cmd))
          logoutput("Saving quality profile plots to quality_Profile_R*.pdf")

          lapply(nameslist, function(x) ggsave(filename=file.path(outdir, paste0("qualityProfile_", x, ".pdf")), plot = pqp[[x]], width=8, height=qplotheight, units="in"))

      })
  }


# Filter and trim ---------------------------------------------------------
  trimlist <- lapply(readslist, function(y) sapply(y, FUN = function(x) gsub(paste0(".", tools::file_ext(x)), "_trim.fastq.gz", x), USE.NAMES = FALSE))

  if (data_type == "PE") {
    cmd <- paste0("out <- filterAndTrim(fwd=file.path(datadir,readslist$R1), filt=file.path(filt.dir,trimlist$R1),rev=file.path(datadir,readslist$R2), filt.rev=file.path(filt.dir,trimlist$R2),  maxEE=", deparse(maxEE),", trimLeft=",deparse(trimLeft),", truncQ=", truncQ, ", truncLen = ",deparse(truncLen) ,", rm.phix=TRUE, compress=TRUE, verbose=TRUE, multithread=nthread, minLen=50)")
  } else {
    cmd <- paste0("out <- filterAndTrim(fwd=file.path(datadir,readslist$R1), filt=file.path(filt.dir,trimlist$R1), maxEE=",maxEE,", trimLeft=",deparse(trimLeft),", truncQ=",truncQ,", truncLen=", deparse(truncLen),", rm.phix=TRUE, compress=TRUE, verbose=TRUE, multithread=nthread, minLen=50)")
  }

  run_cmd(cmd, "dada2::filterAndTrim", bline=1, w2e="No reads passed the filter.")
  print(out)


  ## Check for missing trimmed files
  logoutput("Checking that trimmed files exist.")
  run_cmd('list2env(checktrimfiles(A, filt.dir, trimlist), envir = environment())', "dada2::filterAndTrim")


# Learn errors ------------------------------------------------------------
  cmd <- paste0("err <- lapply(trimlist, function(x) learnErrors(x, multithread=nthread, nbases=", dparams$nbases, ",randomize=TRUE))")
  run_cmd(cmd, "dada2::learnErrors", bline=1)

  run_cmd("pe <- lapply(err, function(x) plotErrors(x, nominalQ=TRUE))", "dada2::plotErrors", bline=1)

  lapply(nameslist, function(x)   ggsave(filename=file.path(outdir, paste0("errorRate_", x, ".pdf")), pe[[x]]))

#  run_cmd("lapply(nameslist, function(x) paste(c(x, dada2:::checkConvergence(err[[x]])), collapse =' ') )", "dada2:::checkConvergence", bline=1)

  trimlist <- lapply(trimlist, function(x) {
    names(x) <- A$SampleID
    return(x)
  })

  saveRDS(err, file = file.path(interm.dir, "err.rds"))
  err <- readRDS(file=file.path(interm.dir, "err.rds"))


#  Dereplicate reads, run dada and merge reads ----------------------------
## loop through each sample individually following: https://benjjneb.github.io/dada2/bigdata.html

  sampleVariants <- foreach (sample=A$SampleID, .inorder=TRUE, .combine='comb', .multicombine=T, .init=list(sv=list(), denoisedF=list(), denoisedR=list(), merged=list())) %do% {

    run_cmd("derep <- lapply(trimlist, function(x) derepFastq(x[sample], verbose=TRUE))", "dada2::derepFastq", bline = 1)

    run_cmd("dd <- sapply(nameslist, function(x) dada(derep[[x]], err=err[[x]], multithread=nthread, verbose=F), USE.NAMES=TRUE, simplify=FALSE)", "dada2::dada")

    ## print information about dada sequence variants
    cat(sapply(nameslist, function(x) { paste(x, capture.output(show(dd[[x]]))[2], sep=": ", collapse='\n') }, USE.NAMES = F))

    if (data_type == "PE") {
      cmd <- paste0("mergePairs(dd$R1, derep$R1, dd$R2, derep$R2, verbose=TRUE, minOverlap=",dparams$minOverlap,", trimOverhang=", trimOverhang,", maxMismatch=", maxMismatch, ", justConcatenate=", justConcatenate, ")")
      sv <- run_cmd(cmd, "dada2::mergePairs", bline=1)
      list(sv=sv, denoisedF=getN(dd$R1), denoisedR=getN(dd$R2), merged=getN(sv))
    } else {
      list(sv=dd$R1, denoisedF=getN(dd$R1), denoisedR=NA, merged=NA)
    }
  }

  names(sampleVariants$sv) <- A$SampleID
  run_cmd( "seqtab <- makeSequenceTable(sampleVariants$sv)", "dada2::makeSequenceTable", bline=1)
  saveRDS(seqtab, file.path(interm.dir, "seqtab.rds"))
  sampleVariants$sv <- NULL


  logoutput(paste0("Removing sequences of length less than ", dparams$min_seq_length, "bp"))
  run_cmd("seqlengths <- nchar(colnames(seqtab))")
  cmd <- paste0("seqtab <- seqtab[,which(seqlengths >=", dparams$min_seq_length ,"), drop=F]")
  run_cmd(cmd)
  cmd <- paste0('saveRDS(seqtab, file.path(interm.dir,"seqtab_min', dparams$min_seq_length, '.rds"))')
  run_cmd(cmd)
  sampleVariants[[ paste0("filter", dparams$min_seq_length)]] <- rowSums(seqtab)



#  Remove chimera if option is chosen -------------------------------------
  if (chimera) {
    run_cmd("seqtabnochimera <- removeBimeraDenovo(seqtab, verbose=TRUE, multithread=nthread)", "dada2::removeBimeraDenovo", bline=1)

    saveRDS(seqtabnochimera, file.path(interm.dir, "seqtab_nochimera.rds"))
    postchimreads <- 100*sum(seqtabnochimera)/sum(seqtab)
    writeLines(paste0("% Reads remaining after chimera removal: ", postchimreads))
    run_cmd("seqtab <- seqtabnochimera")
    sampleVariants$nochim <- rowSums(seqtabnochimera)
  }

  logoutput("Track Reads", 1)
  print(do.call(cbind, lapply(sampleVariants, unlist)))


  ## Write sequence fasta file
  run_cmd('rep_seq_names <- make_seq_names(seqtab)')
  seqs <- dimnames(seqtab)[[2]]
  names(seqs) <- replace_names(seqs, rep_seq_names)
  cmd <- paste0('writeFasta(seqs, file=file.path(outdir, "', dparams$outputfasta,'"))')
  run_cmd(cmd)


# Taxonomic assignment ----------------------------------------------------
  logoutput(paste('Taxonomic assignment with', taxmethod), 1)
  taxtab <- 'taxa'
  if (tolower(taxmethod) == 'idtaxa') {
    taxa <- decipher_assign(refdb = refdb, seqtab = seqtab, nthread=nthread)
    ## OTU table will only go up to genus

  } else {
  ## rdp with assignTaxonomy
    cmd <- paste0("taxa <- assignTaxonomy(seqtab, refdb, multithread=nthread, minBoot=",minBoot,", tryRC=TRUE, verbose=TRUE)")
    run_cmd(cmd, bline=1)

    if (!is.null(refdb_species)) {
      ## OTU table will have species info
      taxtab <- 'taxa.species'

      if (!justConcatenate | data_type == "SE") {
        run_cmd("taxa.species <- addSpecies(taxa, refdb_species, verbose=TRUE, tryRC=TRUE, n=4000)", "dada2::addSpecies", bline=1)
        gc(verbose=T)
      } else {
        ## potentially add this https://github.com/benjjneb/dada2/issues/529
        logoutput("Species assigment only with R1 because justConcatenate is TRUE.", 1)
        ## mergePairs justConcatenate adds 10Ns between reads https://www.bioconductor.org/packages/release/bioc/manuals/dada2/man/dada2.pdf#Rfn.mergePairs.1
        run_cmd('sq.R1 <- str_match(getSequences(seqtab), "(.+?)N{10}")[,2]')
        run_cmd('gen.spc <- assignSpecies(sq.R1, refdb_species, tryRC=T, verbose=TRUE, n=4000)',  "dada2::assignSpecies")
        ## To duplicate addSpecies functionality
        run_cmd('same.genus <- gen.spc[,"Genus"] == taxa[,"Genus"]', "dada2::assignSpecies")
        run_cmd('same.genus[is.na(same.genus)] <- FALSE')
        run_cmd('taxa.species <- cbind(taxa, Species = rep(NA, nrow(taxa)))')
        run_cmd('taxa.species[same.genus, "Species"] <- gen.spc[same.genus, "Species"]')
      }
    }
  }


# Write output files ------------------------------------------------------
  ## Write biom file
  run_cmd('otu_tab <- seqtab; colnames(otu_tab) <- replace_names(colnames(otu_tab), rep_seq_names)', bline=1)
  run_cmd(paste0('row.names(', taxtab, ') <- replace_names(row.names(', taxtab, '), rep_seq_names)'))
  cmd <- paste0('write_biom(dada2biom(otu_tab,' ,taxtab,', metadata = A), file.path(outdir, "',dparams$biomfile,'"))')
  run_cmd(cmd)

  ## Write OTU table
  cmd <- paste0('dada2text(otu_tab, ', taxtab,', file.path(outdir, "', dparams$otutable, '"))')
  run_cmd(cmd)

  ## Write taxonomy table
  cmd <-  paste0('dada2taxonomy(', taxtab,', file.path(outdir, "', dparams$taxtable, '"))')
  run_cmd(cmd)

  ## run garbage collection just in case
  rm(list=setdiff(ls(), lsf.str()))
  gc(verbose=T)

  ## Return 0 on success
  return(as.integer(0))
}

#'
#' wrapper for dada2compute
#'
#' @name trycomputewrapper
#'
#' @rdname dada2compute
#'
#' @description \code{trycomputewrapper} is a wrapper for \code{dada2compute} to be used with
#' rpy2 to run the pipeline in python.  It sets up the global R options and the output to the
#' logfile before running \code{dada2compute}.
#'
#' @param datadir Directory containing FASTQ files.
#' @param outdir Output directory.
#' @param mapfile Mapping filename - needs to be tab separated and have, at minimum,
#'   SampleID, ForwardFastqFile, and ReverseFastqFile columns.
#' @param logfilename Log file name (full path).
#' @param ... parameters to pass to \code{dada2compute}
#'
#' @return \code{trycomputewrapper} returns 0 if \code{dada2compute} succeeds and stops on error.
#' In rpy2, it will throw \code{rpy2.rinterface.RRuntimeError}.
#'
#' @source [computation.R](../R/computation.R#L513)
#'
#' @export
#'
trycomputewrapper <- function(datadir, outdir, mapfile, logfilename="logfile.txt", ...) {
  oldopts <- options()
  options(stringsAsFactors = FALSE, scipen = 999, warn=1, showErrorCalls=TRUE, showNCalls=500, digits.secs = 3, echo=F)

  # dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  logfile <- file(logfilename, open = "at")
  sink(file = logfile, type="output")
  sink(file = logfile, type= "message")

  on.exit(oldopts)
  on.exit(closeAllConnections(), add = TRUE)
  retvalue <- dada2compute(datadir, outdir, mapfile, ...)

  return(retvalue)
}
