% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_iqtree.R
\name{run_iqtree}
\alias{run_iqtree}
\title{IQ-TREE analysis.}
\usage{
run_iqtree(
  command,
  args,
  prefix,
  outdir,
  retain,
  coding = FALSE,
  stdout = "",
  stderr = ""
)
}
\arguments{
\item{command}{a character string with the system command. If not supplied, either "iqtree2" or "iqtree.exe" is used,
depending on the platform. Also the path to directory with the appropriate binary file can be specified
instead of the command itself.}

\item{args}{a character vector in the form accepted by \code{args} argument of \link[base:system2]{base::system2} or a data frame with arguments
to \code{command} (argument names in the 1st column, values in the 2nd column) or a name of a tab-delimited file with the arguments.}

\item{prefix}{a character string with prefix for all output files. It has priority over the specification
through \code{args} argument. If not specified, whether through \code{prefix} or \code{args}, the input file name
without extension is used.}

\item{outdir}{a character string with the path to output directory. It has priority over the path included in prefix.
If not specified in any way, it is assumed to be the same as the working directory.}

\item{retain}{regular expressions with identifiers of outputs to be retained. If specified, the other outputs are deleted.}

\item{coding}{logical, if \code{TRUE} a proper partition file is created and used in the substitution model selection.}

\item{stdout, stderr}{where output to ‘stdout’ or ‘stderr’ should be sent, as in \link[base:system2]{base::system2}.}
}
\description{
Runs IQ-TREE analysis specified by the arguments, including post-processing of the output files.
}
\details{
The minimum \code{args} contains single row with "-s" and the name of the input file. The input directory is
either retrieved fron the input name (if it contains a path) or assummed to be the working directory. The argument
names must be supplied with the leading dashes ("-" or "--", what is appropriate). For the reference to IQ-TREE commands,
consult ... http://www.iqtree.org/doc/Command-Reference ...
}
