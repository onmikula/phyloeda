#' Finding SNPs.
#' 
#' @description
#' Finds variable positions (a.k.a. single nucleotide polymorphisms or SNPs) in sequence alignements.
#'
#' @param loci a list of locus-specific sequence alignments of class `"matrix"` or `"DNAbin"`,
#'   possibly an output of [import_seq].
#' @param allele a regular expression defining an allele identifier (allows to distinguish individual IDs).
#' @returns A list of lists with components `"snp"` (numeric, indices of SNP positions),
#'   `"pis"` (logical, whether the SNPs are parsimony-informative)
#'   and `"nalleles"` (the number of distinct haplotypes = alleles in the alignment).
#' @export

find_snps <- function(loci, allele=NULL) {
	count_states <- function(x) sapply(c("A", "C", "G", "T"), function(s) sum(x == s, na.rm=TRUE))
	if (nchar(system.file(package="dplyr")) > 0) {
		count_alleles <- function(x) dplyr::n_distinct(x, na.rm=TRUE)
	} else {
		count_alleles <- function(x) length(na.omit(unique(x)))
	}
	if (is.null(allele)) {
		is_pis <- function(x, ind) sum(count_states(x) > 1) > 1
	} else {
		is_pis <- function(x, ind) sum(colSums(table(ind, x) > 0) > 1) > 1
	}
	
	searching <- function(locus, modify, allele) {
		if (length(locus) == 0) {
			return(list(snp=integer(0), pis=logical(0), nalleles=integer(0)))
		} else {
			locus <- modify(locus)
			locus[!locus %in% c("A", "C", "G", "T")] <- NA
			nalleles <- apply(locus, 2, count_alleles)
			snp <- which(nalleles > 1)
			if (is.null(allele)) {
				ind <- NULL
			} else {
				ind <- sub(allele, "", rownames(locus))
			}
			if (length(snp) == 0) {
				pis <- integer(0)
			} else {
				pis <- sapply(snp, function(i, ind) is_pis(locus[, i], ind), ind=ind)
			}
			return(list(snp=snp, pis=pis, nalleles=nalleles[snp]))
		}
	}
	is_list <- is.list(loci)
	if (isFALSE(is_list)) {
		loci <- list(loci)
	}
	modify <- ifelse(inherits(loci[[1]], "DNAbin"), base::toupper, base::identity)
	snps <- lapply(loci, searching, modify=modify, allele=allele)
	if (isFALSE(is_list)) {
		snps <- snps[[1]]
	}
	return(snps)
}



#' Extracting SNPs.
#' 
#' @description
#' Extracts specified SNPs from locus-specific alignments.
#'
#' @param loci A list of locus-specific sequence alignments of class `"matrix"` or `"DNAbin"`.
#' @param indices A list of SNP indices produced by [find_snps]. If not provided, it is obtained internally.
#' @param pis Logical, whether to retain only parsimony-informative sites. If allele is specified,
#'   it requires the minority variant to be present in at least two individuals.
#' @param nalleles Numeric, the maximum no. of alleles per variable site,
#'   i.e., `nalleles=2` means extracting just biallelic SNPs.
#' @param rm.invariant Logical, whether to remove invariant loci from the list.
#' @param toupper Logical, whether to enforce representation in upper case letters.
#' @param subset A vector defining subset of loci; logical, numeric (indices of loci) or character (names of loci).
#' @param allele A regular expression defining an allele identifier (allows to distinguish individual IDs).
#' @returns A list of locus-specific SNP matrices (subsets of original alignments), each with attributes `"snp"`,
#'   `"pis"` and`"nalleles"` corresponding to `indices` argument.
#' @export

get_snps <- function(loci, indices=NULL, pis=FALSE, nalleles=4, rm.invariant=TRUE, toupper=FALSE, subset=TRUE, allele=NULL) {
	if (isTRUE(subset)) {
		subset <- rep(TRUE, length(loci))
	} else if (is.numeric(subset)) {
		subset <- seq_along(loci) %in% subset
	} else if  (is.character(subset)) {
		subset <- names(loci) %in% subset
	}
	bedtable <- attr(loci, "bedtable")
	if (!is.null(bedtable)) {
		bedtable <- bedtable[subset,]
	}
	loci <- loci[subset]
	if (!is.null(indices)) {
		indices <- indices[subset]
	}
	extraction <- function(x, snp, modify) {
		x <- modify(x[, snp$snp, drop=FALSE])
		attr(x, "snp") <- snp$snp
		attr(x, "pis") <- snp$pis
		attr(x, "nalleles") <- snp$nalleles
		return(x)
	}
	if (is.null(indices)) {
		indices <- find_snps(loci, allele=allele)
	} else {
		indices <- indices
	}
	if (isTRUE(pis)) {
		indices <- lapply(indices, function(x) list(snp=x$snp[x$pis], pis=x$pis[x$pis], nalleles=x$nalleles[x$pis]))	
	}
	if (nalleles < 4) {
		indices <- lapply(indices, function(x, n) {nall <- x$nalleles <= n; list(snp=x$snp[nall], pis=x$pis[nall], nalleles=x$nalleles[nall])}, n=nalleles)	
	}
	if (isTRUE(rm.invariant)) {
		nz <- sapply(lapply(indices, "[[", "snp"), length) > 0
	} else {
		nz <- rep(TRUE, length(loci))
	}
	modify <- ifelse(toupper == TRUE, base::toupper, base::identity)
	snps <- mapply(extraction, loci[nz], indices[nz], MoreArgs=list(modify=modify), SIMPLIFY=FALSE)
	if (!is.null(bedtable)) {
		attr(snps, "bedtable") <- bedtable[nz,,drop=FALSE]
	}
	return(snps)
}


#' Single SNP per locus.
#' 
#' @description
#' Selecting one SNP per locus, randomly within constraints given by the function arguments.
#'
#' @param snps A list of locus-specific SNP matrices, as produced by [get_snps].
#' @param pis Logical, whether to retain only parsimony-informative sites.
#' @param nalleles Numeric, the maximum no. of alleles per variable site.
#' @param minmiss Logical, whether to select preferentially SNPs with the minimum of mossing data.
#' @param drop Logical, whether to drop out SNPs not complying the criteria (default is `TRUE`), 
#'   or whether to give them lower priority (if `FALSE`).
#' @param seed A seed for random sampling, its specification enforces repeatability of the random selection.
#' @returns A list of locus-specific SNP matrices, but with just one site (column) per locus. Invariant loci a discarded.
#' @export

single_snps <- function(snps, pis=FALSE, nalleles=4, minmiss=FALSE, drop=TRUE, seed=sample(1e+8,1)) {	
	bedtable <- attr(snps, "bedtable")
	pos <- lapply(snps, attr, which="snp")
	allsnps <- lapply(lapply(snps, ncol), seq)
	if (isTRUE(pis) & nalleles < 4) {
		selected <- lapply(snps, function(x) which(attr(x, "pis") & attr(x, "nalleles") <= nalleles))
	} else if (isTRUE(pis)) {
		selected <- lapply(snps, function(x) which(attr(x, "pis")))
	} else if (nalleles < 4) {
		selected <- lapply(snps, function(x) which(attr(x, "nalleles") <= nalleles))
	} else {
		selected <- allsnps
	}
	if (!drop) {
		none_selected <- sapply(selected, length) == 0
		if (any(none_selected)) {
			selected <- Map(function(cond, x, y) {x[which(cond)] <- y[which(cond)]; x}, none_selected, selected, allsnps)	
		}	
	}
	nz <- sapply(selected, length) > 0
	selected <- selected[nz]
	pos <- pos[nz]
	snps <- snps[nz]
	if (!is.null(bedtable)) {
		bedtable <- bedtable[nz,]
	}
	if (isTRUE(minmiss)) {
		nucl <- c("A","C","G","T","a","c","g","t")
		miss <- Map(function(x, col) tapply(!x[,col,drop=FALSE] %in% nucl, rep(seq_along(col), each=nrow(x)), sum), snps, selected)
		selected <- Map(function(col, mis) col[which(mis == min(mis))], selected, miss)
	}
	set.seed(seed)
	selected <- lapply(selected, function(x) ifelse(length(x) == 1, x, sample(x, 1)))
	pos <- Map("[[", pos, selected)
	snps <- Map(function(x, col) x[,col,drop=FALSE], snps, selected)
	if (!is.null(bedtable)) {
		bedtable[,2] <- bedtable[,2] + unlist(pos)
		bedtable[,3] <- bedtable[,2] + 1
	}
	attr(snps, "bedtable") <- bedtable
	return(snps)
}



#' Biallelic SNPs.
#' 
#' @description
#' Makes SNPs biallelic by marking the third and fourth alleles as missing data.
#'
#' @param snps A list of locus-specific SNP matrices, as produced by [get_snps].
#' @param toupper Logical, whether to enforce representation in upper case letters.
#' @returns A list of modified locus-specific SNP matrices.
#' @export

as_biallelic <- function(snps, toupper=FALSE) {
	count_states <- function(x) sapply(c("A", "C", "G", "T"), function(s) sum(x == s, na.rm=TRUE))
	if (isTRUE(toupper)) {
		snps <- lapply(snps, toupper)
	}
	nalleles <- lapply(snps, attr, "nalleles")
	multi <- lapply(nalleles, function(x) which(x > 2))
	coords <- cbind(locus=rep(seq_along(snps), sapply(multi, length)), site=unlist(multi))
	if (nrow(coords) > 0) {
		for (i in seq(nrow(coords))) {
			site <- snps[[coords[i,1]]][,coords[i,2]]
			counts <- count_states(site)
			minor <- sample(c("A", "C", "G", "T"), 4)
			minor <- minor[order(counts[minor])][1:2]
			snps[[coords[i,1]]][site %in% minor,coords[i,2]] <- "N"
			attr(snps[[coords[i,1]]], "nalleles")[coords[i,2]] <- 3
		}
	}
	return(snps)
}



#' Counting SNPs.
#' 
#' @description
#' A utility function counting SNPs specified the function arguments.
#'
#' @param indices A list of SNP indices produced by [find_snps].
#' @param pis Logical, whether to retain only parsimony-informative sites.
#' @param nalleles Numeric, the maximum no. of alleles per variable site.
#' @returns A list of locus-specific SNP matrices, but with just one site (column) per locus. Invariant loci a discarded.
#' @export

count_snps <- function(indices, pis=FALSE, nalleles=4) {
	if (length(indices) == 3 & identical(names(indices), c("snp","pis","nalleles"))) {
		indices <- list(indices)
	}
	if (isTRUE(pis)) {
		indices <- lapply(indices, function(x) list(snp=x$snp[x$pis], x$nalleles[x$pis]))	
	}
	if (nalleles < 4) {
		indices <- lapply(indices, function(x, n) x$snp[x$nalleles <= n], n=nalleles)	
	} else {
		indices <- lapply(indices, "[[", "snp")
	}
	counts <- sapply(indices, length)
	return(counts)
}


