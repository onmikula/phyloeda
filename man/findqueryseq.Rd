% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classifyseq.R
\name{findqueryseq}
\alias{findqueryseq}
\title{Query sequences.}
\usage{
findqueryseq(ref, haps)
}
\arguments{
\item{ref}{a character vector with IDs of reference sequences or a single tree (readable by \link{read_tree}),
whose tip labels are used as IDs of reference sequences.}

\item{haps}{either an output of \link{haplotypes} or a data frame with assignment of sequences to haplotypes
(\code{assign} component of the list produced by \link{haplotypes} or a name of tab-delimited file with such data).}
}
\value{
A character vector with IDs of query sequences.
}
\description{
Finds query sequences for phylogenetic placement, i.e., those that are different from reference sequences and from each other.
}
