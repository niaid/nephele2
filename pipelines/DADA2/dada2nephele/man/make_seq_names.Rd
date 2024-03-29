% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ioutils.R
\name{make_seq_names}
\alias{make_seq_names}
\alias{replace_names}
\title{ASV sequence names}
\source{
[ioutils.R](../R/ioutils.R#L18)
}
\usage{
make_seq_names(seqtab, type = "simple")

replace_names(tabnames, seq_names)
}
\arguments{
\item{seqtab}{OTU/DADA2 Sequence table}

\item{type}{one of: "simple" renaming seq1, seq2, ...
or "md5" for hash like in QIIME 2}

\item{tabnames}{character of original names}

\item{seq_names}{named vector with elements being new names and names from
\code{tabnames}, e.g. output of \emph{make_seq_names}}
}
\value{
\emph{make_seq_names} returns named vector where elements are names and the names are the
sequences

\emph{replace_names} returns (unnamed) vector of new names
}
\description{
DADA2 uses the actual ASV sequences as names for sequence and
taxonomy tables which is unwieldy for the output.  We make simpler new names
seq1, seq2, ...

\emph{make_seq_names} makes the new names

\emph{replace_names} takes in old names and returns new ones.
}
