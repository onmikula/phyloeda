#' Concatenate sequences.
#' 
#' @description
#' Concatenates multi-locus sequences.
#' 
#' @param loci a list of alignments in `matrix` or `DNAbin` format.
#' @param toupper logical, whether to export sequences in upper case letters.
#' @param missing a symbol representing missing nucleotides.
#' @returns A concatenated alignment in `matrix` or `DNAbin` format with `part` attribute defining locus boundaries.
#' @export

concatenate <- function(loci, toupper=TRUE, missing="-") {
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
		as.DNAbin <- TRUE
	} else {
		as.DNAbin <- FALSE
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



#' Partition sequences.
#' 
#' @description
#' Partitions multi-locus sequences.
#' 
#' @param x a sequence alignment in `matrix` or `DNAbin` format.
#' @param part either a two column matrix indicating starts and ends of loci in the concatenated aligments
#'  (in .fasta or .nexus files) or a name of partition file formatted according to RAxML or MrBayes convention.
#' @param toupper logical, whether to export sequences in upper case letters.
#' @param missing a symbol representing missing nucleotides.
#' @returns A list of alignments in `matrix` or `DNAbin` format with `part` attribute defining locus boundaries.
#' @export

partition <- function(x, part=NULL, toupper=TRUE) {
	if (is.null(part)) {
		part <- attr(x, "part")
	} else if (is.character(part) & length(part) == 1) {
		part <- read_part(part)
	} else if (is.matrix(part) | is.data.frame(part)) {
		part <- as.matrix(part[,1:2])
	}
	if (inherits(x, "DNAbin")) {
		x <- as.character(x)
		as.DNAbin <- TRUE
	} else {
		as.DNAbin <- FALSE
	}
	if (is.null(part)) {
		part <- matrix(seq(ncol(x)), ncol(x), 2, dimnames=list(colnames(x), c("start", "end")))
	} else {
		colnames(part) <- c("start", "end")
	}
	loci <- setNames(vector("list", nrow(part)), rownames(part))
	for (i in seq_along(loci)) {
		loci[[i]] <- x[,part[i,1]:part[i,2],drop=FALSE]
	}
	if (isTRUE(toupper)) {
		loci <- lapply(loci, toupper)
	}
	if (isTRUE(as.DNAbin)) {
		loci <- lapply(loci, ape::as.DNAbin)
	}
	attr(loci, "part") <- part
	return(loci)
}

