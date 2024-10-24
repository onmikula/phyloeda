% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_str_output.R
\name{custom_str_palette}
\alias{custom_str_palette}
\title{Custom STRUCTURE palette.}
\usage{
custom_str_palette(out, info, palette)
}
\arguments{
\item{out}{data frame produced by \link{read_str_output}.}

\item{info}{a data frame listing individual labels (1st column) and population labels (2nd column).}

\item{palette}{a named vector of colors, the names are labels of populations from \code{info}.}
}
\description{
Prepares a custom palette for \link{plot_str_output}.
}
