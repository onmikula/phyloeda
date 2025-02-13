% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc_traces.R
\name{combine_runs}
\alias{combine_runs}
\title{Combining posterior samples.}
\usage{
combine_runs(x)
}
\arguments{
\item{x}{a list of \code{multiPhylo} objects or a data frames.}
}
\value{
An object of class \code{multiPhylo} or a data frame containing pooled posterior samples.
}
\description{
Combines posterior samples obtained in independent runs of the same analysis.
}
