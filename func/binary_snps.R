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
		snps <- apply(snps, 2, function(n) match(n, rep_len(names(rev(sort(table(n)))),2))) - 1
	}
	mode(snps) <- "numeric"
	rownames(snps) <- rows
	if (isTRUE(center)) {
		snps <- snps - 0.5
	}
	if (isTRUE(onerowperind) & all(grepl(allele, rownames(snps)))) {
		snps <- do.call(rbind, split(snps, gsub(allele, "", rownames(snps))))
		snps <- do.call(cbind, by(t(snps), rep(seq(ncol(snps) / 2), each=2), colSums))
	}
#	if (isTRUE(onerowperind) & all(!grepl(allele, rownames(x)))) {}
#	if (isFALSE(onerowperind) & all(!grepl(allele, rownames(x)))) {}
	return(snps)
}


### FUNCTION
# impute_binary_snps
### ARGUMENTS
# b: the output of make_binary_snps
# d: matrix genotype dissimilarities
# k: no. of nearest neighbors used for imputation of missing data

impute_binary_snps <- function(b, d=NULL, k=3) {
	if (is.null(d)) {
		d <- dissim.snps(b)
	} else if (inherits(d, "dist")) {
		d <- as.matrix(d)
	}
	nn <- apply(d, 1, order, simplify=FALSE)
	nn <- Map(setdiff, nn, as.list(seq_along(nn)))
	for (i in seq(nrow(b))) {
		nas <- which(is.na(b[i,]))
		for (j in nas) {
			a <- na.omit(b[nn[[i]],j])
			if (length(a) > 0) {
				b[i,j] <- mean(a[1:min(c(k, length(a)))])
			}			
		}
	}
	return(b)
}



### FUNCTION
# subset_binary_snps
### ARGUMENTS
# b: the output of make_binary_snps
# rows:
# cols:

subset_binary_snps <- function(b, rows=NULL, cols=NULL) {
	rsubs <- function(b) apply(!is.na(b), 1, all)
	csubs <- function(b) apply(!is.na(b), 2, all)
	if (is.null(rows) & is.null(cols)) {
		nr <- nrow(b)
		nc <- ncol(b)
		b <- b[rsubs(b), csubs(b)]
		while (nrow(b) != nr & ncol(b) != nc) {
			nr <- nrow(b)
			nc <- ncol(b)
			b <- b[rsubs(b), csubs(b)]	
		}
	} else if (is.null(rows)) {
		b <- b[rsubs(b), cols]
	} else if (is.null(cols)) {
		b <- b[rows, csubs(b)]		
	} else {
		b <- b[rows, cols]
	}
	return(b)
}



### FUNCTION
# center_binary_snps
### ARGUMENTS
# b: the output of make_binary_snps

center_binary_snps <- function(b, onerowperind=TRUE) {
	return(b - ifelse(onerowperind, 1, 0.5))
}



### FUNCTION
# dissim.snps
### ARGUMENTS
# b: the output of make_binary_snps

dissim.snps <- function(b) {
	hamming <- function(x, y) {
		d <- abs(x - y)
		return(sum(d, na.rm=TRUE) / sum(!is.na(d)))
	}
	b <- b / diff(range(b, na.rm=TRUE))
	dst <- matrix(0, nrow(b), nrow(b), dimnames=list(rownames(b), rownames(b)))
	for (i in 1:(nrow(b) - 1)) {
		for (j in (i+1):nrow(b)) {
			dst[i,j] <- dst[j,i] <- hamming(b[i,], b[j,])
		}
	}
	return(dst)
}

