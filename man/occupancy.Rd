% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loci_utils.R
\name{occupancy}
\alias{occupancy}
\title{Occupancy.}
\usage{
occupancy(loci)
}
\arguments{
\item{loci}{a list of locus-specific sequence alignments of class \code{"matrix"} or \code{"DNAbin"}.}
}
\value{
A numeric vector of locus-specific occupancies.
}
\description{
Calculates occupancy, i.e., proportion of individuals for which a locus was successfully sequenced.
}
