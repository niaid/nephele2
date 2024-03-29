% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ioutils.R
\name{dada2output}
\alias{dada2output}
\alias{dada2biom}
\alias{dada2text}
\alias{dada2taxonomy}
\title{dada2output}
\source{
[ioutils.R](../R/ioutils.R#L60)
}
\usage{
dada2biom(otu, tax, metadata = NULL)

dada2text(otu, tax, filename)

dada2taxonomy(tax, filename)
}
\arguments{
\item{otu}{OTU table.}

\item{tax}{Taxonomy table}

\item{metadata}{Metadata table (Optional).}

\item{filename}{Output filename}
}
\value{
dada2biom returns a biom object.
}
\description{
Utilities to convert dada2 sequence tables to output files.

dada2biom makes valid biom object.

dada2text writes tab separated OTU table to filename.

dada2output writes tab separated taxonomy file in format
suitable for importing into qiime2
}
