library(devtools)
library(rmarkdown)
library(Rd2md)
options(warn=1)
devtools::load_all()


# Functions and constants -------------------------------------------

addlinenumber <- function(filename) {
  cfile <- readLines(filename)
  ll <- grep("@source", cfile)
  if (length(ll) == 0) return(NULL)
  newl <- sapply(ll, function(x) gsub("\\.R\\)|\\.R#L\\d+\\)", paste0(".R#L", x, ")"), cfile[x], perl=T))
  cfile[ll] <- newl
  writeLines(cfile, con=filename)
  return(ll)
}

myoutputoptions <- list(pandoc_args=c("--template", file.path(find.package("rmarkdown"), "rmarkdown/templates/github_document/resources/default.md"), "--atx-headers"))

remove.file <- function(x) {
  suppressWarnings(file.remove(x))
}


clean_pandoc2_highlight_tags = function(x) {
  # x = gsub('(</a></code></pre>)</div>', '\\1', x)
  x = gsub('<div class="sourceCode" id="cb\\d+">', '', x)
  x = gsub('</div></li>', '</li>', x, perl=T)
  x = gsub('</pre></div>', '</pre>', x, perl=T)
  # x = gsub('<a class="sourceLine"[^>]+>(.*)</a>', '\\1', x)
  # x = gsub('(<code class="sourceCode .+?">)<span id=".+?"><a href=".+?"></a>', '\\1', x)
  # x = gsub('</span>(</code></pre>)', '\\1', x)
  x = gsub('<span id="cb.*?"><a href="#cb.*?"></a>', '',x, perl=TRUE)
  x = gsub('</span>$', '', x, perl=T)
  x
}

# Build package -----------------------------------------------------------

## Add line number to source lines in roxygen comments
Rfiles <- list.files("./R", full.names = T)
sapply(Rfiles, addlinenumber)

## document
document(roclets=c('rd', 'collate', 'namespace'))

## Imported packages - can check DESCRIPTION
ns <- scan("NAMESPACE", sep="\n", what = character())
importedpackages <- unique(stringr::str_match(ns, "import.*\\((.*?)[\\,\\)]")[,2])

## Description
desc::desc_normalize()

## install
devtools::install(args = c("--preclean", "--no-multiarch", "--with-keep.source", "--no-staged-install"), dependencies = FALSE, upgrade=FALSE)
# devtools::load_all()

stop()
# Make dada2_config.py
options(useFancyQuotes=FALSE)
# writeparams <- getOption("dparams")
# writeparams$trimLeft <- deparse(writeparams$trimLeft)
# writeLines(c('#!/usr/bin/env python3', '## Automatically generated.  Do Not Edit.  Default parameters.  See nephele2/pipelines/DADA2/dada2nephele/R/computation.R .onLoad', unlist(lapply(names(writeparams), function(x) paste0(toupper(x), '=', sQuote(writeparams[[x]]))))), con = '../dada2_config.py')
#

# Package Documentation ---------------------------------------------------


## library function specification =====

# mdfile <- "Reference_Manual_dada2nephele.md"
# Rmdfile <- gsub(".md", ".Rmd", mdfile)
#
# ReferenceManual(outdir = file.path(getwd(), "doc"), man_file = Rmdfile,  title.level = 2, run.examples = TRUE, sepexported = TRUE, toc.matter = NULL, code.headings = F, topic.section.heading = F)
# render_manual_github(file.path("doc", Rmdfile), outdir = file.path(getwd(), "doc"), toc=T, toc_depth = 3, knitr_opts_chunk = list(echo=T, tidy.opts = list(width.cutoff=80)), toplinks = F, nocodelinks=T, author='Poorani Subramanian')
# file.remove(file.path("doc", Rmdfile))
#
#



# User docs ---------------------------------------------------------------

userRmd <- "doc/user_doc.Rmd"
mdfile <- "doc/user_doc.md"
render(userRmd, output_format = "github_document", output_options=myoutputoptions)

Sys.setenv(PATH="/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/opt/X11/bin:/Users/subramanianp4/Applications/RStudio.app/Contents/MacOS/postback:/Users/subramanianp4/Library/Python/3.9/bin:/Users/subramanianp4/Library/Nodejs/bin")
render("doc/user_doc.Rmd", output_format = "html_document", output_file = "dada2_pipeline.html", output_options=list(pandoc_args = c("--ascii", "-F", "panflute")), clean=T)
htmlfile <- readLines("doc/dada2_pipeline.html")
htmlfile <- sapply(htmlfile, clean_pandoc2_highlight_tags)

writeLines(htmlfile, "doc/dada2_pipeline.html")
jsb <- system2(command='which', args='js-beautify', stdout=T)
if (is.null(attr(jsb, "status"))) system2(jsb, args=c('--type', 'html', '-b', 'collapse', '-n', '-r', '-f', "doc/dada2_pipeline.html"))
if (file.exists("doc/dada2_pipeline_files")) unlink("doc/dada2_pipeline_files", recursive = T)


###     template: ../../../../resources/Rdocs/user_doc_template.html
