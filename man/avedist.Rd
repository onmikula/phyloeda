% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/avedist.R
\name{avedist}
\alias{avedist}
\title{Average genetic distances.}
\usage{
avedist(
  loci,
  info = NULL,
  seqmap = NULL,
  allele = NULL,
  model = "raw",
  diag = TRUE
)
}
\arguments{
\item{loci}{The list of locus-specific alignments or a single such alignment (possibly
with concatenated multi-locus data). The row names must be interpretable as individual labels
(possibly with help of \code{seqmap} or \code{allele} arguments).}

\item{info}{A matrix or data frame (1st column individuals, 2nd column species)
or a name of a tab-delimited file with such data.}

\item{seqmap}{A data frame mapping names of sequences (1st column) to those of individuals (2nd column).}

\item{allele}{An allele identifier, regular expression distinguishing sequences from the same individual.}

\item{model}{The nucleotide substitution model indicated as in \code{\link[ape:dist.dna]{ape::dist.dna()}}.}

\item{diag}{Logical, whether to calculate average intraspecific distances, default is \code{TRUE}.}
}
\value{
A square symmetric matrix giving mean genetic distance between species (off-diagonal)
and within them (diagonal).
}
\description{
Calculates average genetic distance between the specified species (populations, lineages) and within them.
}
