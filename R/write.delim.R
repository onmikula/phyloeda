#' Writing tab-delimited file.
#' 
#' @description
#' A wrapper for [utils::write.table] writing a tab-delimited text file. The counterpart of [utils::read.delim].
#'
#' @param x a matrix or data frame to be written.
#' @param file character string, name of the output file.
#' @param row.names logical, whether to use row names (default is `FALSE`).
#' @param col.names logical, whether to use column names (default is `TRUE`).
#' @param quote logical or numeric, whether to enclose field contents into quotes (or in which columns).
#' @param sep the field separator string.
#' @export

write.delim <- function(x, file, row.names=FALSE, col.names=TRUE, quote=FALSE) {
	utils::write.table(x, file, row.names=row.names, col.names=col.names, quote=quote, sep="\t")
}
