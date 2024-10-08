% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cartoon_tree.R
\name{cartoon_tree}
\alias{cartoon_tree}
\title{Cartoon tree.}
\usage{
cartoon_tree(
  phy,
  info,
  color,
  alpha = 0.5,
  lwd = 2,
  cex = 1,
  show.sp.label = TRUE,
  type = "phylogram",
  use.edge.length = TRUE,
  direction = "rightwards",
  device,
  file,
  width = 7,
  height = 7,
  mai = 0.02,
  xylim = 4,
  return = FALSE
)
}
\arguments{
\item{phy}{an object of class \code{phylo}.}

\item{color}{named vector of colors, the names correspond to clade labels.
Optionally, an element called \code{"bg"} can be included to define color of unclassified branches.}

\item{alpha}{numeric, opacity of triangle filling.}

\item{lwd}{width of the lines.}

\item{cex}{size of species labels.}

\item{show.sp.label}{logical, whether to show species labels.}

\item{type}{character, tree type, either \code{"phylogram"} (default) or \code{"cladogram"}.}

\item{use.edge.length}{logical, whether to use edge lengths.}

\item{direction}{a character string specifyying direction of the tree.}

\item{device}{Either \code{"quartz"}, \code{"x11"} or \code{"pdf"} to indicate the type of graphical device.
If \code{NULL}, the objects are plotted into the current device.}

\item{file}{a name of pdf file if \code{device == "pdf"}.}

\item{width}{width of graphical device in inches.}

\item{height}{height of graphical device in inches.}

\item{mai}{size of outer margins in inches, recycled if necessary.}

\item{xylim}{percents of node coordinate ranges to be added to x and y axes at their extremes.}

\item{return}{logical, whether to return coordinates of species labels.}

\item{cartoon}{a data frame (or matrix) with tip labels (1st column) and clade assignment (2nd column).
The classification need not to be comprehensive.}
}
\description{
Plots phylogenetic tree while collapsing defined clades into triangles.
}
