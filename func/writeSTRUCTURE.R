writeSTRUCTURE <- function(loci, file=NULL, return=FALSE) {
	alleles <- vector("list", length(loci))
	for (i in seq_along(loci)) {
		bp <- matrix(toupper(loci[[i]]) %in% c("A","C","G","T"), nrow(loci[[i]]), ncol(loci[[i]]))
		sq <- loci[[i]][,apply(bp, 2, all), drop=FALSE]
		sq <- apply(sq, 1, paste, collapse="")
		sq <- sq[order(rownames(loci[[i]]))]
		ind <- gsub("_[[:digit:]]+$", "", rownames(loci[[i]]))
		al <- setNames(match(sq, unique(sq)), ind)
		alleles[[i]] <- matrix(al, length(al)/2, 2, byrow=TRUE, dimnames=list(unique(ind), NULL))
	}
	rows <- sort(unique(unlist(lapply(alleles, names))))
	structure <- rep(list(matrix(-9, length(rows), 2, dimnames=list(rows, NULL))), length(alleles))
	for (i in seq_along(structure)) {
		structure[[i]][match(rownames(al), rows),] <- alleles[[i]]
	}
	structure <- do.call(cbind, structure)
	if (!is.null(file)) {
		write.table(structure, file=file, col.names=FALSE, row.names=FALSE, quote=FALSE)
	}
	if (isTRUE(return)) {
		return(structure)
	}
}


