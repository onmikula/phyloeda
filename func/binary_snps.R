### FUNCTION
# make_binary_snps
### DESCRIPTION
# converts biallelic snps to their binary representation
### ARGUMENTS
# snps: SNPs in 'matrix' of 'DNAbin' format of R package 'ape'
# center: whether to represent genotypes as centered, i.e. with 0 zero for heterozygote and -1, 1 for homozygotes
# allele: allele identifier in the rownames of snps
# onerowperind: whether to put individual genotypes on a single row 
# format: data format - either 'snps' (default) or 'structure', diifers in treatment of ambiguous and missing data

make_binary_snps <- function(snps, center=FALSE, scale=FALSE, onerowperind=TRUE, allele="_[[:digit:]]$", format=c("snps","structure")) {
	if (is.list(snps)) {
		bp <- cumsum(sapply(snps, ncol))
		part <- cbind(c(1, bp[-length(bp)] + 1), unname(bp))
		rows <- sort(unique(unlist(lapply(snps, rownames))))
		concatenated <- matrix(NA, length(rows), max(bp), dimnames=list(rows, NULL))
		for (i in seq_along(snps)) {
			concatenated[rownames(snps[[i]]), part[i,1]:part[i,2]] <- snps[[i]]
		}
		snps <- concatenated
	}	
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
		ord <- match(rownames(b), rownames(d))
		d <- d[ord, ord]
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
# rows, cols: indices of rows and columns to be retained

subset_binary_snps <- function(b, rows=NULL, cols=NULL) {
	if (is.null(rows)) {
		rows <- rep(TRUE, nrow(b))
	}
	if (is.null(cols)) {
		cols <- rep(TRUE, ncol(b))
	}
	return(b[rows, cols])
}



### FUNCTION
# center_binary_snps
### ARGUMENTS
# b: the output of make_binary_snps

center_binary_snps <- function(b, onerowperind=TRUE) {
	return(b - ifelse(onerowperind, 1, 0.5))
}


### FUNCTION
# scale_binary_snps
# scales genotypes by the expected rate of genetic drift
### ARGUMENTS
# b: the output of make_binary_snps

scale_binary_snps <- function(b, center=TRUE, onerowperind=TRUE) {
	multi <- ifelse(isTRUE(onerowperind), 2, 1)
	p <- colSums(b, na.rm=TRUE) / (multi * colSums(!is.na(b)))
	b <- base::scale(b, center=TRUE, scale=FALSE)
	bcenter <- colMeans(b, na.rm=TRUE)
	bscale <- sqrt(p * (1 - p))
	for (j in seq(ncol(b))) {
		b[,j] <- (b[,j] - bcenter[j]) / bscale[j]
	}
	if (isFALSE(center)) {
		b <- b + rep(1, nrow(b)) %*% t(bcenter)
	}
	return(b)
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



### FUNCTION
# write.binary
### ARGUMENTS
# b: the output of make_binary_snps
# file: file name

write.binary <- function(b, file) {
	write.table(b, file, row.names=TRUE, col.names=FALSE, quote=FALSE, sep="\t")
}


### FUNCTION
# read.binary
### ARGUMENTS
# file: file name

read.binary <- function(file) {
	as.matrix(read.table(file, header=FALSE, row.names=1, sep="\t"))
}

