#' Concatenate sequences.
#' 
#' @description
#' Concatenates multi-locus sequences.
#' 
#' @param loci a list of alignments in `matrix` or `DNAbin` format.
#' @param toupper logical, whether to export sequences in upper case letters.
#' @param as.DNAbin logical, whether to return the alignment in `DNAbin` format.
#' @param missing a symbol representing missing nucleotides.
#' @returns A concatenated alignment in `matrix` or `DNAbin` format with `part` attribute defining locus boundaries.
#' @export

concatenate <- function(loci, toupper=TRUE, as.DNAbin=FALSE, missing="-") {
	bedtable <- attr(loci, "bedtable")
	bp <- cumsum(sapply(loci, ncol))
	part <- attr(loci, "part")
	if (is.null(part)) {
		part <- cbind(c(1, bp[-length(bp)] + 1), unname(bp))
		rownames(part) <- names(bp)
	}
	if (!is.null(attr(loci[[1]], "snp"))) {
		snps <- lapply(loci, function(x) cbind(snp=attr(x, "snp"), pis=as.character(attr(x, "pis")), nalleles=attr(x, "nalleles")))
		snps <- cbind(locus=rep(names(loci), sapply(snps, nrow)), do.call(rbind, snps))
	} else {
		snps <- NULL
	}
	if (inherits(loci[[1]], "DNAbin")) {
		loci <- lapply(loci, as.character)
	}
	rows <- sort(unique(unlist(lapply(loci, rownames))))
	concatenated <- matrix(missing, length(rows), max(bp), dimnames=list(rows, NULL))
	for (i in seq_along(loci)) {
		concatenated[rownames(loci[[i]]), part[i,1]:part[i,2]] <- loci[[i]]
	}
	if (isTRUE(toupper)) {
		concatenated <- toupper(concatenated)
	}
	if (isTRUE(as.DNAbin)) {
		concatenated <- ape::as.DNAbin(concatenated)
	}
	attr(concatenated, "snps") <- snps
	attr(concatenated, "bedtable") <- bedtable
	attr(concatenated, "part") <- part
	return(concatenated)
}
