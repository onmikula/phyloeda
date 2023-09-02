# Function:
#	write_diem_input
# Arguments:
#	snps: a matrix or a list of matrices with SNP genotypes, rows are individual alleles, columns are SNP sites
#	file: name of output file
#	return: whether to return the output as a matrix
#	allele: allele identifier in the row names

write_diem_input <- function(snps, file=NULL, return=FALSE, allele="_[012]$", missing="N") {
	code_states <- function(x) {
		n <- na.omit(unique(c("A", "C", "G", "T", x)))
		cx <- sapply(n, function(s) sum(x == s, na.rm=TRUE))
		cx <- c(sort(cx[1:4], decreasing=TRUE), cx[-seq(4)])
		cz <- min(c(3, which(cx == 0)))
		if (cz > 2) {
			cx <- replace(cx, 1, 0)
			cx <- replace(cx, 2, 1)
			cx <- replace(cx, 3:length(cx), NA)
		} else {
			cx <- replace(cx, 1:length(cx), NA)
		}
		return(cx)
	}
	if (is.list(snps)) {
		bp <- cumsum(sapply(snps, ncol))
		part <- attr(snps, "part")
		if (is.null(part)) {
			part <- cbind(c(1, bp[-length(bp)] + 1), unname(bp))
			rownames(part) <- names(bp)
		}
		snps <- lapply(snps, as.matrix)
		rows <- sort(unique(unlist(lapply(snps, rownames))))
		concatenated <- matrix(missing, length(rows), max(bp), dimnames=list(rows, NULL))
		for (i in seq_along(snps)) {
			concatenated[rownames(snps[[i]]), part[i,1]:part[i,2]] <- snps[[i]]
		}
		snps <- concatenated
	}
	snps <- toupper(snps)
	snps[!snps %in% c("A", "C", "G", "T") ] <- NA
	states <- apply(snps, 2, code_states, simplify=FALSE)
	for (i in seq(ncol(snps))) {		
		snps[,i] <- states[[i]][snps[,i]]
	}
	snps <- snps[,colSums(is.na(snps)) < nrow(snps),drop=FALSE]
	mode(snps) <- "numeric"
	indivs <- gsub(allele, "", rownames(snps))
	counts <- do.call(cbind, by(snps, indivs, colSums, na.rm=TRUE))[,unique(indivs)]
	unknown <- do.call(cbind, by(is.na(snps), indivs, colSums, na.rm=TRUE))[,unique(indivs)] > 0
	counts[unknown] <- "U"
	mode(counts) <- "character"
	if (!is.null(file)) {
		writeLines(paste0(apply(cbind("S", counts), 1, paste, collapse="")), file)
	}
	if (isTRUE(return)) {
		return(counts)
	}
}
