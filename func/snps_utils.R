### FUNCTION
# find_snps
### ARGUMENTS
# loci: list of locus-specific alignments in 'matrix' or 'DNAbin' format of R package 'ape'
# nalleles: the maximum no. of alleles per site (i.e., nalleles=2 extracts just biallelic SNPs)

find_snps <- function(loci, nalleles=4) {
	searching <- function(loci, nalleles) {
		if (length(loci) == 0) {
			return(c(var=0, pis=0))
		} else {
			loci <- toupper(as.matrix(loci))
			loci[! loci %in% c("A", "C", "G", "T") ] <- NA
			ndist <- apply(loci, 2, dplyr::n_distinct, na.rm=TRUE)
			vs <- which(ndist > 1 & ndist <= nalleles)
			is <- which(sapply(vs, function(i) sum(table(loci[, i]) > 1)) > 1)
			return(list(var=vs, pis=vs[is]))
		}
	}
	if (is.list(loci)) {
		snps <- lapply(loci, searching, nalleles=nalleles)
	} else {
		snps <- searching(loci, nalleles=nalleles)
	}
	return(snps)
}


### FUNCTION
# count_snps 
### ARGUMENTS
# snps: list of SNP indices produced by find_snps
# type: either 'pis' (parsimony informative sites) or 'var' (all variable sites) or both (default)

count_snps <- function(snps, type=c("var","pis")) {
	if (length(snps) == 2 & identical(names(snps), c("var","pis"))) {
		snps <- list(snps)
	}
	counts <- sapply(lapply(snps, function(x) lapply(x, length)), "[[", type[1])
	return(counts)
}


### FUNCTION
# get_snps
### DESCRIPTION
# extracts SNPs from a list of locus-specific alignments in 'matrix' or 'DNAbin' format of R package 'ape'
### ARGUMENTS
# loci: list of locus-specific alignments in 'matrix' of 'DNAbin' format of R package 'ape'
# nalleles: the maximum no. of alleles per site (i.e., nalleles=2 extracts just biallelic SNPs)
# type: either 'pis' (parsimony informative sites) or 'var' (all variable sites)
# single: retain just a single randomly chosen SNP per locus
# as.matrix: whether to produce alignments in 'matrix' format

get_snps <- function(loci, nalleles=4, type="pis", single=FALSE, as.matrix=TRUE) {
	snps_pos <- lapply(find_snps(loci, nalleles= nalleles), "[[", type)
	nz <- sapply(snps_pos, length) > 0
	loci <- loci[nz]
	snps_pos <- snps_pos[nz]
	if (isTRUE(single)) {
		snps_pos <- lapply(snps_pos, function(x) ifelse(length(x) == 1, x, sample(x, 1)))
	}
	snps <- Map(function(x, cols) x[,cols, drop=FALSE], loci, snps_pos)
	snps <- lapply(lapply(snps, as.matrix), toupper)	
	if (as.matrix == TRUE) {
		rows <- sort(unique(unlist(lapply(snps, rownames))))
		snps <- do.call(cbind, lapply(snps, function(x, r) x[match(r, rownames(x)),], r=rows))
		snps[is.na(snps)] <- "N"
	}
	return(snps)
}


### FUNCTION
# bpcoverage 
### ARGUMENTS
# x: a locus-specific alignement

bpcoverage <- function(x) {
	x[!toupper(x) %in% c("A", "C", "G", "T") ] <- NA
	return(colSums(!is.na(x)) / nrow(x))
}



### FUNCTION
# occupancy 
### ARGUMENTS
# loci: list of locus-specific alignments

occupancy <- function(loci) {
	seqnam <- sort(unique(unname(unlist(lapply(loci, rownames)))))
	occup <- sapply(loci, nrow) / length(seqnam)
	return(occup)
}



### FUNCTION
# locilengths 
### ARGUMENTS
# loci: list of locus-specific alignements

locilengths <- function(loci) {
	put_nas <- function(x) { x[!x %in% c("A", "C", "G", "T") ] <- NA; return(x) }
	omit_nas <- function(x) x[, !apply(is.na(x), 2, all), drop=FALSE]
	if (is.matrix(loci)) {
		loci <- list(loci)
	}
	loci <- lapply(lapply(loci, as.matrix), toupper)
	loci <- lapply(lapply(loci, put_nas), omit_nas)
	lengths <- sapply(loci, ncol)
	return(lengths)
}






