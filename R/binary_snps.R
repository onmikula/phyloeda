#' Binary coding of biallelic SNPs.
#' 
#' @description
#' Converts biallelic SNPs to their binary representation.
#' 
#' @param snps single nucleotide polymorphisms in `matrix` or `DNAbin` format.
#' @param allele regular expression, allele identifier in the rownames of snps.
#' @param counts logical, whether to represent individual genotype in a single row by counts of '1' alleles.
#' @param format data format, either `"snps"` (default) or `"structure"`, diifers in treatment of ambiguous and missing data.
#' @param center logical, whether to represent genotypes as centered, i.e. with 0 zero for heterozygote and -1, 1 for homozygotes.
#' @param scale logical, whether to scale genotypes by the expected rate of genetic drift. If `center=FALSE` and `scale=TRUE`,
#'   the data are also centered so the column means are zero. By default, both `center` and `scale` are `FALSE`.
#' @returns A matrix with genotypes represented by `0` and `1` alleles or their counts (if `counts=TRUE`).
#' @export

make_binary_snps <- function(snps, allele="_[[:digit:]]$", counts=TRUE, format=c("snps","structure"), center=FALSE, scale=FALSE) {
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
	if (isTRUE(counts) & all(grepl(allele, rownames(snps)))) {
		snps <- do.call(rbind, split(snps, gsub(allele, "", rownames(snps))))
		snps <- do.call(cbind, by(t(snps), rep(seq(ncol(snps) / 2), each=2), colSums))
	}
	if (isTRUE(scale)) {
		snps <- scale_binary_snps(snps, center=!center, counts=counts)
	}
	
	return(snps)
}


#' Imputation of biallelic SNPs.
#' 
#' @description
#' Imputes missing biallelic snps.
#' 
#' @param b matrix, the output of [make_binary_snps].
#' @param d optional, matrix of pairwise genotype dissimilarities.
#' @param k numeric, the number of nearest neighbor genotypes used for imputation of missing data.
#' @returns A modified matrix with binary coded genotypes.
#' @export

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



#' Subset of SNPs.
#' 
#' @description
#' Subsets binary coded SNPs.
#' 
#' @param b matrix, the output of [make_binary_snps].
#' @param rows, cols indices of rows and columns to be retained.
#' @returns A matrix with subset of binary coded genotypes.
#' @export

subset_binary_snps <- function(b, rows=NULL, cols=NULL) {
	if (is.null(rows)) {
		rows <- rep(TRUE, nrow(b))
	}
	if (is.null(cols)) {
		cols <- rep(TRUE, ncol(b))
	}
	return(b[rows, cols])
}



#' Center SNP genotypes.
#' 
#' @description
#' Centers binary coded SNP genotypes, so in the 'counts' representation, heterozygote takes value of 0 and homozygotes -1 and 1.
#' 
#' @param b matrix, the output of [make_binary_snps].
#' @param counts logical, whether individual genotype is represented in a single row by counts of '1' alleles.
#' @returns A matrix with centered binary coded genotypes.
#' @export

center_binary_snps <- function(b, counts=TRUE) {
	return(b - ifelse(counts, 1, 0.5))
}


#' Scale SNP genotypes.
#' 
#' @description
#' Scales and possibly also centers binary coded SNP genotypes. The centering is different from [center_binary_snps],
#'   as it just makes all column means to be zero.
#' 
#' @param b matrix, the output of [make_binary_snps].
#' @param center logical, whether to apply centering (default is `TRUE`).
#' @param counts logical, whether individual genotype is represented in a single row by counts of '1' alleles.
#' @returns A matrix with centered binary coded genotypes.
#' @export

scale_binary_snps <- function(b, center=TRUE, counts=TRUE) {
	multi <- ifelse(isTRUE(counts), 2, 1)
	p <- colSums(b, na.rm=TRUE) / (multi * colSums(!is.na(b)))
	b <- base::scale(b, center=center, scale=FALSE)
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


#' Dissimilarity of SNP genotypes.
#' 
#' @description
#' Calculates matrix of pairwise genotype dissimilarities. The dissimilarity is just number of different alleles,
#'   averaged over alleles and omitting data missing in either of the two compared genotypes.
#' 
#' @param b matrix, the output of [make_binary_snps].
#' @param allele regular expression, allele identifier in the rownames of snps.
#' @returns A matrix of pairwise dissimilarities.
#' @export

dissim.snps <- function(b, allele="_[[:digit:]]$") {
	hamming <- function(x, y) {
		d <- abs(x - y)
		return(sum(d, na.rm=TRUE) / sum(!is.na(d)))
	}
	if (all(grepl(allele, rownames(b)))) {
		indnam <- sub(allele, "", rownames(b))
		b <- do.call(rbind, by(b, indnam, colSums, na.rm=TRUE))
		b <- b[order(match(rownames(b), indnam)),]
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



#' Writing file with SNP genotypes.
#' 
#' @description
#' Writes a tab-delimited file with binary coded SNP genotypes.
#' 
#' @param b matrix, the output of [make_binary_snps].
#' @param file character string, the file name.
#' @export

### FUNCTION
# write.binary
### ARGUMENTS
# b: the output of make_binary_snps
# file: file name

write_binary <- function(b, file) {
	write.table(b, file, row.names=TRUE, col.names=FALSE, quote=FALSE, sep="\t")
}


#' Reading file with SNP genotypes.
#' 
#' @description
#' Read a tab-delimited file with binary coded SNP genotypes.
#' 
#' @param file character string, the file name.
#' @returns A matrix with binary coded SNP genotypes.
#' @export

read_binary <- function(file) {
	as.matrix(read.table(file, header=FALSE, row.names=1, sep="\t"))
}

