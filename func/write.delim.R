write.delim <- function(x, file, row.names=FALSE, col.names=TRUE, quote=FALSE) {
	write.table(x, file, row.names=row.names, col.names=col.names, quote=quote, sep="\t")
}
