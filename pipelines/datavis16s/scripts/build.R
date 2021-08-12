library(devtools)
library(rmarkdown)
library(Rd2md)

Sys.setenv(PATH=paste(Sys.getenv("PATH"), file.path(system2("python3", c("-m", "site","--user-base"), stdout=T), "bin"), sep=":"))

# Build package ----------------------------------------------------------------------------------
document(roclets=c('rd', 'collate', 'namespace'))
devtools::install(args = c("--preclean", "--no-multiarch", "--with-keep.source"), dependencies = F, upgrade=F)



## Imported packages - can check DESCRIPTION
ns <- scan("NAMESPACE", sep="\n", what = character())
importedpackages <- unique(stringr::str_match(ns, "import.*\\((.*?)[\\,\\)]")[,2])


## update description
desc::desc_set(Date=format(Sys.time(), format="%F"), normalize=TRUE)
deptable <- desc::desc_get_deps()
# deptable$version <- apply(deptable, 1, function(x) { if (x[2] == "R") return(x[3]); paste("==", packageVersion(x[2])) })
desc::desc_set_deps(deptable, normalize = TRUE)


# Documentation --------------------------------

myoutputoptions <- list(pandoc_args=c("--template", file.path(find.package("rmarkdown"), "rmarkdown/templates/github_document/resources/default.md"), "--atx-headers", "--columns=10000"))


remove.file <- function(x) {
  suppressWarnings(file.remove(x))
}

# clean_pandoc2_highlight_tags = function(x) {
#   x = gsub('(</a></code></pre>)</div>', '\\1', x)
#   x = gsub('<div class="sourceCode"[^>]+>(<pre)', '\\1', x)
#   x = gsub('<a class="sourceLine"[^>]+>(.*)</a>', '\\1', x)
#   x
# }

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
#  x = gsub('<div id="[\\w-]+" class="mt-5 mb-5">', '<div class="mt-5 mb-5">', x, perl=T)
  x
}



## library function specification ============

mdfile <- "Reference_Manual_datavis16s.md"
Rmdfile <- gsub(".md", ".Rmd", mdfile)
## CRAN Rd2md
# ReferenceManual(outdir = file.path(getwd(), "doc"), front.matter = yaml)
## My Rd2md https://github.com/pooranis/Rd2md
ReferenceManual(outdir = file.path(getwd(), "doc"), man_file = Rmdfile,  title.level = 2, run.examples = FALSE, sepexported = TRUE, toc.matter = NULL, code.headings = F, topic.section.heading = F)
render_manual_github(file.path("doc", Rmdfile), outdir = file.path(getwd(), "doc"), toc=T, toc_depth = 3, knitr_opts_chunk = list(echo=T, eval=F, tidy.opts = list(width.cutoff=80)), nocodelinks=T)
file.remove(file.path("doc", Rmdfile))

stop()
# # User docs ---------------------------------------------------------------
#
# userRmd <- "doc/user_doc.Rmd"
# mdfile <- "doc/user_doc.md"
# oops <- myoutputoptions
# oops$pandoc_args <- c(oops$pandoc_args, "-M", "title=User docs")
# render(userRmd, output_format = "md_document", output_options = oops)

# htmlfilename <- "doc/datavis16s_pipeline.html"
# render("doc/user_doc.Rmd", output_format = "html_document", output_file = basename(htmlfilename), output_options=list(pandoc_args = c("--ascii", "-F", "panflute")))
# htmlfile <- readLines(htmlfilename)
# htmlfile <- sapply(htmlfile, clean_pandoc2_highlight_tags)
# writeLines(htmlfile, htmlfilename)
# jsb <- system2(command='which', args='js-beautify', stdout=T)
# if (is.null(attr(jsb, "status"))) system2(jsb, args=c('--type', 'html', '-n', '-r', '-f', htmlfilename))






