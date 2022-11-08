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
	snps <- lapply(lapply(lapply(snps, as.matrix), as.character), toupper)	
	if (as.matrix == TRUE) {
		rows <- sort(unique(unlist(lapply(snps, rownames))))
		snps <- do.call(cbind, lapply(snps, function(x, r) x[match(r, rownames(x)),], r=rows))
		snps[is.na(snps)] <- "N"
	}
	return(snps)
}


### FUNCTION
# make_binary_snps
### DESCRIPTION
# converts biallelic snps to their binary representation
### ARGUMENTS
# snps: SNPs in 'matrix' of 'DNAbin' format of R package 'ape'
# format: data format - either 'snps' (default) or 'structure', diifers in treatment of ambiguous and missing data
# center: whether to represent genotypes as centered, i.e. with 0 zero for heterozygote and -1, 1 for homozygotes
# allele: allele identifier in the rownames of snps
# onerowperind: whether to put individual genotypes on a single row 

make_binary_snps <- function(snps, format=c("snps","structure"), center=FALSE, onerowperind=TRUE, allele="_[[:digit:]]$") {
	snps <- toupper(as.matrix(snps))
	rows <- rownames(snps)
	if (format[1] == "structure") {
		snps[snps %in% c("N","-","-9")] <- NA
		loci <- rep(1:(ncol(snps)/2), each=2)
		for (i in 1:(ncol(snps)/2)) {
			xi <- as.vector(snps[,loci == i])
			snps[,loci == i] <- match(xi, sort(unique(xi))) - 1
		}
	} else if (format[1] == "snps") {
		snps[!snps %in% c("A", "C", "G", "T") ] <- NA
		snps <- apply(snps, 2, function(n) match(n, sort(unique(n)))) - 1
	}
	mode(snps) <- "numeric"
	rownames(snps) <- rows
	if (isTRUE(onerowperind) & all(grepl(allele, rownames(snps)))) {
		snps <- do.call(rbind, split(snps, gsub(allele, "", rownames(snps))))
	}
# TBD!	if (isFALSE(onerowperind) & all(!grepl(allele, rownames(x)))) {}
	if (isTRUE(center)) {
		snps <- 2 * (snps - 0.5)
		snps[is.na(snps)] <- 0 
	}
	return(snps)
}



### FUNCTION
# count_binary_snps
### DESCRIPTION
# counts the SNPs site by site (treating every site as independent)
### ARGUMENTS
# b: binary representation of SNPs produced by make_binary_snps
### DETAILS
# currently it assumes diploidy and onerowperind=TRUE 

count_binary_snps <- function(b) {
	center <- any(b == -1)
	if (center) {
		b <- b / 2 + 0.5
		b[b == 0.5] <- NA
	}
	nloci <- ncol(b) / 2
	b <- do.call(cbind, by(t(b), rep(seq(nloci), each=2), colSums))
	if (center) {
		b <- 2 * (b - 1)
		b[is.na(b)] <- 0 
	}
	return(b)
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






