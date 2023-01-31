### FUNCTION
# find_snps
### ARGUMENTS
# loci: list of locus-specific alignments in 'matrix' or 'DNAbin' format of R package 'ape'

find_snps <- function(loci) {
	count_states <- function(x) sapply(c("A", "C", "G", "T"), function(s) sum(x == s, na.rm=TRUE))
	searching <- function(locus, modify) {
		if (length(locus) == 0) {
			return(list(snp=integer(0), pis=logical(0), nalleles=integer(0)))
		} else {
			locus <- modify(locus)
			locus[!locus %in% c("A", "C", "G", "T")] <- NA
			nalleles <- apply(locus, 2, dplyr::n_distinct, na.rm=TRUE)
#			nalleles <- apply(locus, 2, function(x) length(na.omit(unique(x))))
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
### ARGUMENTSz# snps: list of SNP indices produced by find_snps
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

get_snps <- function(loci, snps=NULL, pis=FALSE, nalleles=4, toupper=FALSE) {
	extraction <- function(x, snp, modify) {
		x <- modify(x[, snp$snp, drop=FALSE])
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
	nz <- sapply(lapply(snps_pos, "[[", "snp"), length) > 0
	modify <- ifelse(toupper == TRUE, base::toupper, base::identity)
	snps <- mapply(extraction, loci[nz], snps_pos[nz], MoreArgs=list(modify=modify), SIMPLIFY=FALSE)
	return(snps)
}


### FUNCTION
# single_snps
### DESCRIPTION
# extracts a single SNP per locus
### ARGUMENTS
# snps: list of SNPs produced by get_snps
# pis: whether to retain only parsimony-informative sites
# nalleles: retain sites with the maximum no. of alleles <= nalleles
# drop: whether to drop out SNPs not complying the criteria or to give them lower priority, default is TRUE (drop them out)
# minmiss: whether to retain sites with the minimum of missing data (after applying 'pis' and 'nalleles' criteria)
# seed: seed for random sampling

single_snps <- function(snps, pis=FALSE, nalleles=4, minmiss=FALSE, drop=TRUE, seed=sample(1e+8,1)) {	
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
	snps <- snps[nz]
	if (isTRUE(minmiss)) {
		nucl <- c("A","C","G","T","a","c","g","t")
		miss <- Map(function(x, col) tapply(!x[,col,drop=FALSE] %in% nucl, rep(seq_along(col), each=nrow(x)), sum), snps, selected)
		selected <- Map(function(col, mis) col[which(mis == min(mis))], selected, miss)
	}
	set.seed(seed)
	selected <- lapply(selected, function(x) ifelse(length(x) == 1, x, sample(x, 1)))
	snps <- Map(function(x, col) x[,col,drop=FALSE], snps, selected)
	return(snps)
}



### FUNCTION
# as_biallelic
### DESCRIPTION
# makes SNPs biallelic by marking the third and fourth allels as missing data
### ARGUMENTS
# snps: list of SNPs produced by get_snps
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
# loci: list of locus-specific alignments

locilengths <- function(loci) {
	sapply(loci, ncol)
}



### FUNCTION
# clean_data 
### ARGUMENTS
# loci: list of locus-specific alignments or a single such alignment
# remove: "empty" (default) means removing all empty rows and columns, whereas "rows" and "columns" mean removing rows or columns with any missing data (=retaining only complete ones)

clean_data <- function(loci, remove="empty", numeric=FALSE) {
	cleaning <- function(x, remove, numeric) {
		if (isTRUE(numeric)) {
			notna <- !is.na(x)
		} else {
			notna <- matrix(as.character(x) %in% c("A", "C", "G", "T", "a", "c", "g", "t"), nrow(x), ncol(x))	
		}
		if (remove != "empty") {
			if (remove == "rows") {
				rows <- rowSums(notna) == ncol(x)
				notna <- notna[rows,]
				cleaned <- x[rows,]
				attrib <- c(attributes(cleaned), attributes(x))
				attributes(cleaned) <- attrib[match(names(attributes(x)), names(attrib))]
			}
			if (remove == "columns") {
				cols <- colSums(notna) == nrow(x)
				notna <- notna[,cols]
				cleaned <- x[,cols]
				attrib <- c(attributes(cleaned), attributes(x))
				snpsattr <- intersect(c("pis","nalleles"), names(attrib))
				for (a in snpsattr) {
					attrib[[a]] <- attrib[[a]][cols]
				}
				attributes(cleaned) <- attrib[match(names(attributes(x)), names(attrib))]
			}
		} else if (remove == "empty") {
			rows <- rowSums(notna) > 0
			cols <- colSums(notna) > 0
			cleaned <- x[rows, cols, drop=FALSE]
			attrib <- c(attributes(cleaned), attributes(x))
			snpsattr <- intersect(c("pis","nalleles"), names(attrib))
			for (a in snpsattr) {
					attrib[[a]] <- attrib[[a]][cols]
			}
			attributes(cleaned) <- attrib[match(names(attributes(x)), names(attrib))]
		}
		return(cleaned)	
	}
	is_list <- is.list(loci)
	if (isFALSE(is_list)) {
		loci <- list(loci)
	}
	numeric <- isTRUE(numeric) | mode(loci[[1]]) == "numeric"
	loci <- lapply(loci, cleaning, remove=remove, numeric=numeric)
	if (isFALSE(is_list)) {
		loci <- loci[[1]]
	}
	return(loci)
}



### FUNCTION
# subset_loci 
### ARGUMENTS
# loci: list of locus-specific alignments or a concatenated matrix of such alignments (possibly with 'part' attribute)
# n: number of loci to be selected
# seed: seed for random sampling
# subset: a vector with names or numeric indices of loci or logical values

subset_loci <- function(loci, n, seed=sample(1e+8,1), subset=NULL) {
	is_matrix <- inherits(loci, "matrix")
	nloci <- ifelse(is_matrix, ncol(loci), length(loci))
	if (!is.null(subset)) {
		if (is.logical(subset) & length(subset) != nloci) {
			error("the length of a logical subset must be equal to the number of loci")
		}
	} else {
		if (nloci <= n) {
			warning("n is not smaller than the number of loci, all loci are retained")
			subset <- seq(nloci)
		} else {
			set.seed(seed)
			subset <- sort(sample(nloci, n))
		}
	}
	if (is_matrix) {
		part <- attr(loci, "part")
		if (!is.null(part)) {
			part <- part[subset,]
			cs <- cumsum(part[,2] - part[,1] + 1)
			part[,] <- c(1, cs[-length(cs)] + 1, cs)
		}
		output <- loci[,subset]
		attr(output, "part") <- part
	} else {
		output <- loci[subset]
	}
	return(output)
}



### FUNCTION
# subset_indiv 
### ARGUMENTS
# loci: list of locus-specific alignments or a single such alignment
# n: number of loci to be selected

subset_indiv <- function(loci, subset, allele=NULL) {
	is_matrix <- inherits(loci, "matrix")
	nloci <- ifelse(is_matrix, ncol(loci), length(loci))
	if (is.factor(subset)) {
		subset <- as.character(subset)
	}
	if (is_matrix) {
		part <- attr(loci, "part")
		if (is.character(subset)) {
			nams <- rownames(loci)
			if (!is.null(allele)) {
				nams <- gsub(allele, "", nams)
				subset <- gsub(allele, "", subset)
			}
			loci <- loci[nams %in% subset,]	
		} else {
			loci <- loci[subset,]
		}
		attr(loci, "part") <- part
	} else {
		stopifnot(is.character(subset), error("for the list of loci, 'subset' must by in the form of character vector"))
		if (!is.null(allele)) {
			subset <- gsub(allele, "", subset)
		} else {
			allele <- ""
		}
		for (i in seq_along(loci)) {
			nams <- gsub(allele, "", rownames(loci[[i]]))
			origattrib <- attributes(loci[[i]])
			loci[[i]] <- loci[[i]][nams %in% subset,]
			allattrib <- c(attributes(loci[[i]]), origattrib)
			attributes(loci[[i]]) <- allattrib[match(names(origattrib), names(allattrib))]
		}
	}
	return(loci)	
}

