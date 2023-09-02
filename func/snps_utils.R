
### FUNCTION
# find_snps
### DESCRIPTION
# finds variable positions (SNPs) in the alignement and annotates them
### ARGUMENTS
# loci: list of locus-specific alignments in 'matrix' or 'DNAbin' format of R package 'ape'
### VALUE
# a (list of) list(s) with components 'snp' (numeric, indices of SNP positions), 'pis' (logical, whether are the SNPs parsimony-informative) and 'nalleles' (the number of distinct haplotypes = alleles in the alignment)

find_snps <- function(loci) {
	count_states <- function(x) sapply(c("A", "C", "G", "T"), function(s) sum(x == s, na.rm=TRUE))
	if (nchar(system.file(package="dplyr")) > 0) {
		count_alleles <- function(x) dplyr::n_distinct(x, na.rm=TRUE)
	} else {
		count_alleles <- function(x) length(na.omit(unique(x)))
	}	
	searching <- function(locus, modify) {
		if (length(locus) == 0) {
			return(list(snp=integer(0), pis=logical(0), nalleles=integer(0)))
		} else {
			locus <- modify(locus)
			locus[!locus %in% c("A", "C", "G", "T")] <- NA
			nalleles <- apply(locus, 2, count_alleles)
			snp <- which(nalleles > 1)
			pis <- sapply(snp, function(i) sum(count_states(locus[, i]) > 1)) > 1
			return(list(snp=snp, pis=pis, nalleles=nalleles[snp]))
		}
	}
	is_list <- is.list(loci)
	if (isFALSE(is_list)) {
		loci <- list(loci)
	}
	modify <- ifelse(inherits(loci[[1]], "DNAbin"), base::toupper, base::identity)
	snps <- lapply(loci, searching, modify=modify)
	if (isFALSE(is_list)) {
		snps <- snps[[1]]
	}
	return(snps)
}


### FUNCTION
# count_snps 
### DESCRIPTION
# 
### ARGUMENTS
# snps: list of SNP indices produced by find_snps
# pis: whether to count only parsimony-informative sites
# nalleles: whether to count only parsimony-informative sites

count_snps <- function(snps, pis=FALSE, nalleles=4) {
	if (length(snps) == 3 & identical(names(snps), c("snp","pis","nalleles"))) {
		snps <- list(snps)
	}
	if (isTRUE(pis)) {
		snps <- lapply(snps, function(x) list(snp=x$snp[x$pis], x$nalleles[x$pis]))	
	}
	if (nalleles < 4) {
		snps <- lapply(snps, function(x, n) x$snp[x$nalleles <= n], n=nalleles)	
	} else {
		snps <- lapply(snps, "[[", "snp")
	}
	counts <- sapply(snps, length)
	return(counts)
}



### FUNCTION
# get_snps
### DESCRIPTION
# extracts SNPs from a list of locus-specific alignments in 'matrix' or 'DNAbin' format of R package 'ape'
### ARGUMENTS
# loci: list of locus-specific alignments in 'matrix' of 'DNAbin' format of R package 'ape'
# snps: list of SNP indices produced by find_snps
# pis: whether to retain only parsimony-informative sites
# nalleles: the maximum no. of alleles per site (i.e., nalleles=2 extracts just biallelic SNPs)
# toupper: whether to enforce representation in upper case letters
# subset: subset of loci (logical, numeric - indicesof locior character - names of loci)

get_snps <- function(loci, snps=NULL, pis=FALSE, nalleles=4, toupper=FALSE, rm.invariant=TRUE, subset=TRUE) {
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
	if (!is.null(snps)) {
		snps <- snps[subset]
	}
	extraction <- function(x, snp, modify) {
		x <- modify(x[, snp$snp, drop=FALSE])
		attr(x, "snp") <- snp$snp
		attr(x, "pis") <- snp$pis
		attr(x, "nalleles") <- snp$nalleles
		return(x)
	}
	if (is.null(snps)) {
		snps_pos <- find_snps(loci)
	} else {
		snps_pos <- snps
	}
	if (isTRUE(pis)) {
		snps_pos <- lapply(snps_pos, function(x) list(snp=x$snp[x$pis], pis=x$pis[x$pis], nalleles=x$nalleles[x$pis]))	
	}
	if (nalleles < 4) {
		snps_pos <- lapply(snps_pos, function(x, n) {nall <- x$nalleles <= n; list(snp=x$snp[nall], pis=x$pis[nall], nalleles=x$nalleles[nall])}, n=nalleles)	
	}
	if (isTRUE(rm.invariant)) {
		nz <- sapply(lapply(snps_pos, "[[", "snp"), length) > 0
	} else {
		nz <- rep(TRUE, length(loci))
	}
	modify <- ifelse(toupper == TRUE, base::toupper, base::identity)
	snps <- mapply(extraction, loci[nz], snps_pos[nz], MoreArgs=list(modify=modify), SIMPLIFY=FALSE)
	if (!is.null(bedtable)) {
		attr(snps, "bedtable") <- bedtable[nz,,drop=FALSE]
	}
	return(snps)
}


### FUNCTION
# filter_snps
### DESCRIPTION
# filter SNPs fulfilling specified criteria
### ARGUMENTS
# snps: list of SNPs produced by get_snps
# pis: whether to retain only parsimony-informative sites
# nalleles: retain sites with the maximum no. of alleles <= nalleles

filter_snps <- function(snps, pis=FALSE, nalleles=4, rm.invariant=TRUE) {
	subset_snps <- function(x, col) {
		snp <- attr(x, "snp")
		pis <- attr(x, "pis")
		nalleles <- attr(x, "nalleles")
		x <- x[,col,drop=FALSE]
		attr(x, "snp") <- snp[col]
		attr(x, "pis") <- pis[col]
		attr(x, "nalleles") <- nalleles[col]
		return(x)
	}
	bedtable <- attr(snps, "bedtable")
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

	if (isTRUE(rm.invariant)) {
		nz <- sapply(selected, length) > 0
		selected <- selected[nz]
		snps <- snps[nz]
		if (!is.null(bedtable)) {
			bedtable <- bedtable[nz,,drop=FALSE]
		}
	}

	snps <- Map(subset_snps, snps, selected)
	attr(snps, "bedtable") <- bedtable
	return(snps)
}



### FUNCTION
# single_snps
### DESCRIPTION
# chooses a single SNP per locus
### ARGUMENTS
# snps: list of SNPs produced by get_snps
# pis: whether to retain only parsimony-informative sites
# nalleles: retain sites with the maximum no. of alleles <= nalleles
# minmiss: whether to retain sites with the minimum of missing data (after applying 'pis' and 'nalleles' criteria)
# drop: whether to drop out SNPs not complying the criteria (TRUE, the default) or to give them lower priority
# seed: seed for random sampling

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



### FUNCTION
# as_biallelic
### DESCRIPTION
# makes SNPs biallelic by marking the third and fourth allels as missing data
### ARGUMENTS
# snps: list of SNPs produced by get_snps
# toupper: whether to enforce representation in upper case letters

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



### FUNCTION
# snp_bedtable
### DESCRIPTION
# 
### ARGUMENTS
# snps: list of SNPs produced by get_snps

snp_bedtable <- function(snps) {
	bedtable <- attr(snps, "bedtable")
	snptable <- lapply(snps, attr, which="snp")
	bedtable <- bedtable[rep(seq(nrow(bedtable)), sapply(snptable, length)),]
	snptable <- unname(unlist(snptable))
	bedtable[,2] <- bedtable[,2] + snptable
	bedtable[,3] <- bedtable[,2] + 1	
	bedtable[,4] <- paste0(bedtable[,4], "_pos", snptable)
	rownames(bedtable) <- NULL
	return(bedtable)
} 
