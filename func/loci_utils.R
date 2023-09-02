### FUNCTION
# occupancy: sequencing success across individuals
### ARGUMENTS
# loci: list of locus-specific alignments
### VALUE
# vector with locus-specific occupancies

occupancy <- function(loci) {
	seqnam <- sort(unique(unname(unlist(lapply(loci, rownames)))))
	occup <- sapply(loci, nrow) / length(seqnam)
	return(occup)
}



### FUNCTION
# locilengths: length of loci (including gaps and parts with missing data only)
### ARGUMENTS
# loci: list of locus-specific alignments
# rm.gaps: whether to exclude gaps from calculating of lengths
# outgroup: identifier(s) of outgroup(s) not to be taken into account when detecting gaps
### VALUE
# vector with lengths of loci ()

locilengths <- function(loci, rm.gaps=FALSE, outgroup=NULL) {
	rm_gaps <- function(dna, outgroups) {
		gaps <- toupper(dna[!rownames(dna) %in% outgroups,,drop=FALSE])
		gaps[!gaps %in% c("A","C","G","T")] <- NA
		gaps <- apply(is.na(gaps), 2, all) 
		return(dna[,!gaps,drop=FALSE])	
	}
	if (isTRUE(rm.gaps)) {
		seqnam <- sort(unique(unlist(lapply(loci, rownames))))
		outgroups <- sort(unlist(lapply(outgroup, grep, x=seqnam, value=TRUE)))
		loci <- lapply(loci, rm_gaps, outgroups=outgroups)
	}
	return(sapply(loci, ncol))
}



### FUNCTION
# bpcoverage
### ARGUMENTS
# x: a locus-specific alignement
### VALUE
# vector with proportions of non-missing data along the locus

bpcoverage <- function(x) {
	x[!x %in% c("A", "C", "G", "T", "a", "c", "g", "t") ] <- NA
	bpc <- colSums(!is.na(x)) / nrow(x)
	return(bpc)
}



### FUNCTION
# clean_data
### DESCRIPTION
# removes missing data from the alignment(s)
### ARGUMENTS
# loci: list of locus-specific alignments or a single such alignment
# remove: "empty" (default) means removing all empty rows and columns, whereas "rows" and "columns" mean removing rows or columns with any missing data (=retaining only complete ones)
# numeric: whether are the data numeric (e.g. binary representation of biallelic SNPs)

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
				notna <- notna[rows,,drop=FALSE]
				cleaned <- x[rows,,drop=FALSE]
				attrib <- c(attributes(cleaned), attributes(x))
				attributes(cleaned) <- attrib[match(names(attributes(x)), names(attrib))]
			}
			if (remove == "columns") {
				cols <- colSums(notna) == nrow(x)
				notna <- notna[,cols,drop=FALSE]
				cleaned <- x[,cols,drop=FALSE]
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
### DESCRIPTION
# makes a subset of loci
### ARGUMENTS
# loci: list of locus-specific alignments or a concatenated matrix of such alignments with 'part' attribute
# subset: a vector with names or numeric indices of loci or logical values
# n: number of loci to be selected
# seed: seed for random sampling

subset_loci <- function(loci, subset=NULL, n, seed=sample(1e+8,1)) {
	is_matrix <- inherits(loci, "matrix")
	nloci <- ifelse(is_matrix, nrow(attr(loci, "part")), length(loci))
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
	if (is.character(subset)) {
		subset <- names(loci) %in% subset
	}
	if (is_matrix) {
		bedtable <- attr(loci, "bedtable")
		part <- attr(loci, "part")[subset,,drop=FALSE]
		cols <- unlist(Map(seq, part[,1], part[,2]))
		loci <- loci[,cols,drop=FALSE]
		cs <- cumsum(part[,2] - part[,1] + 1)
		part[,] <- c(1, cs[-length(cs)] + 1, cs)
		attr(loci, "part") <- part
		if (!is.null(bedtable)) {
			attr(loci, "bedtable") <- bedtable[subset,,drop=FALSE]
		}
	} else {
		bedtable <- attr(loci, "bedtable")
		loci <- loci[subset]
		attr(loci, "bedtable") <- bedtable[subset,,drop=FALSE]
	}
	return(loci)
}



### FUNCTION
# subset_indiv: 
### DESCRIPTION
# makes a subset of individuals from the data, while preserving attributes of the alignement(s)
### ARGUMENTS
# loci: list of locus-specific alignments or a single such alignment
# subset: individual indices
# allele: allele identifier

subset_indiv <- function(loci, subset, allele="_[012]$") {

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
		if (!is.character(subset)) {
			stop("for the list of loci, 'subset' must by in the form of character vector")
		}
		if (!is.null(allele)) {
			subset <- sub(allele, "", subset)
		} else {
			allele <- ""
		}
		totorigattrib <- attributes(loci)
		for (i in seq_along(loci)) {
			nams <- sub(allele, "", rownames(loci[[i]]))
			origattrib <- attributes(loci[[i]])
			loci[[i]] <- loci[[i]][nams %in% subset,,drop=FALSE]
			allattrib <- c(attributes(loci[[i]]), origattrib)
			attributes(loci[[i]]) <- allattrib[match(names(origattrib), names(allattrib))]
		}
		totallattrib <- c(attributes(loci), totorigattrib)
		attributes(loci) <- totallattrib[match(names(totorigattrib), names(totallattrib))]
	}
	
	return(loci)	
}
