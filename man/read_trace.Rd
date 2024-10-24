% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc_traces.R
\name{read_trace}
\alias{read_trace}
\title{Reading MCMC traces.}
\usage{
read_trace(file, comment.char = "#", row.names = NULL, sep = "", header = TRUE)
}
\arguments{
\item{file}{name of the input file.}

\item{comment.char}{the commenting character.}

\item{row.names, sep, header}{optional arguments of \link[utils:read.table]{utils::read.table},
relevant if the the imput file is in tabular format.}
}
\value{
An object of class \code{multiPhylo} or a data frame, depending on the nature of the analysis.
}
\description{
Reads posterior sample of trees or parameters obtained in MCMC-based phylogenetic inference.
}
