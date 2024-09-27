#' Occupancy.
#' 
#' @description
#' Calculates occupancy, i.e., proportion of individuals for which a locus was successfully sequenced.
#'
#' @param loci a list of locus-specific sequence alignments of class `"matrix"` or `"DNAbin"`.
#' @returns A numeric vector of locus-specific occupancies.
#' @export

occupancy <- function(loci) {
	seqnam <- sort(unique(unname(unlist(lapply(loci, rownames)))))
	occup <- sapply(loci, nrow) / length(seqnam)
	return(occup)
}



#' Group occupancy.
#' 
#' @description
#' Calculates group occupancy, i.e., proportion of groups (species, populations, lineages ...)
#' for which a locus was successfully sequenced from at least single individual.
#'
#' @param loci a list of locus-specific sequence alignments of class `"matrix"` or `"DNAbin"`.
#' @param info a data frame or matrix listing individual labels (1st column) and group labels (2nd column).
#' @param allele regular expression, allele identifier in the rownames of snps.
#' @returns A numeric vector of locus-specific group occupancies.
#' @export

groccupancy <- function(loci, info, allele="_[012]$") {
	seqnam <- sort(unique(unname(unlist(lapply(loci, rownames)))))
	groups <- info[match(sub(allele, "", seqnam), info[,1]),2]
	ngroup <- length(unique(group))
	groccup <- unname(sapply(loci, function(x) length(unique(groups[match(rownames(x), seqnam)]))) / ngroup)
	return(groccup)
}



#' Coverage.
#' 
#' @description
#' Calculates coverage, i.e., proportion of loci for which an individual was successfully sequenced.
#'
#' @param loci a list of locus-specific sequence alignments of class `"matrix"` or `"DNAbin"`.
#' @param allele a regular expression defining an allele identifier (allows to distinguish individual IDs).
#' @param indiv an optional character vector giving complete list of individual IDs (useful to get also
#'   individuals with zero coverage represented in the output).
#' @param prop logical, whether to express the coverages as proportions (default) or counts
#' @returns A numeric vector of individual-specific coverages.
#' @export

coverage <- function(loci, allele="_[012]$", indiv=NULL, prop=TRUE) {
	succ <- table(unlist(lapply(loci, function(x) unique(sub(allele, "", rownames(x))))))
	if (isTRUE(prop)) {
		succ <- succ / length(loci)
	}
	succ <- setNames(as.numeric(succ), names(succ))
	if (!is.null(indiv)) {
		miss <- setdiff(indiv, names(succ))
		if (length(miss) > 0) {
			succ <- c(succ, setNames(rep(0, length(miss)), miss))
		}
		succ <- succ[indiv]
	}
	return(succ)
}



#' Group coverage.
#' 
#' @description
#' Calculates group coverage, i.e., proportion of loci for which at least one individual from given groups
#'   (species, populations, lineages ...) was successfully sequenced.
#'
#' @param loci a list of locus-specific sequence alignments of class `"matrix"` or `"DNAbin"`.
#' @param info a data frame or matrix listing individual labels (1st column) and group labels (2nd column).
#' @param allele a regular expression defining an allele identifier (allows to distinguish individual IDs).
#' @returns A numeric vector of group-specific coverages.
#' @export

grcoverage <- function(loci, info, allele="_[012]$") {
	indnam <- sort(unique(sub(allele, "", unlist(lapply(loci, rownames)))))
	groups <- info[match(indnam, info[,1]),2]
	succ <- lapply(loci, function(x) unique(groups[match(sub(allele, "", rownames(x)), indnam)]))
	succ <- table(unlist(succ))
	succ <- succ / length(loci)
	succ <- setNames(as.numeric(succ), names(succ))
	return(succ)
}



#' Locus lengths.
#' 
#' @description
#' Calculates length of loci, possibly including gaps and data-defficient parts.
#'
#' @param loci a list of locus-specific sequence alignments of class `"matrix"` or `"DNAbin"`.
#' @param rm.gaps logical, whether to exclude gaps from calculating of lengths.
#' @param outgroup identifier(s) of outgroup(s) not to be taken into account when detecting gaps.
#' @returns A numeric vector of locus lengths.
#' @export

loclen <- function(loci, rm.gaps=FALSE, outgroup=NULL) {
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
	len <- sapply(loci, ncol)
	return(len)
}



#' Base-pair occupancy.
#' 
#' @description
#' For every base pair it calculates its occupancy, i.e., proportion of sequences with non-missing data.
#'
#' @param locus a locus-specific alignment.
#' @returns A numeric vector with proportions of non-missing data along the locus.
#' @export

bpcoverage <- function(locus) {
	locus[! locus %in% c("A", "C", "G", "T", "a", "c", "g", "t") ] <- NA
	bpc <- colSums(!is.na(locus)) / nrow(locus)
	return(bpc)
}



#' Cleaning of data.
#' 
#' @description
#' Removes missing data for (an) alignment(s).
#'
#' @param loci a list of locus-specific sequence alignments of class `"matrix"` or `"DNAbin"`.
#' @param remove `"empty"` (default) means removing all empty rows and columns, which are completely empty,
#'   whereas `"rows"` and `"columns"` means removing whole rows or columns with any missing data.
#' @param numeric logical, whether are the data numeric (e.g. binary representation of biallelic SNPs).
#' @returns A numeric vector of locus-specific occupancies.
#' @export

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



#' Subsets of loci.
#' 
#' @description
#' Makes a subset of loci.
#'
#' @param loci a list of locus-specific sequence alignments of class `"matrix"` or `"DNAbin"`.
#' @param subset a vector with names or numeric indices of loci or logical values.
#' @param n a number of loci to be randomly selected.
#' @param seed numeric, a seed for random sampling.
#' @returns A numeric vector of locus-specific occupancies.
#' @export

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



#' Subsets of individuals.
#' 
#' @description
#' Makes a subset of individuals acroiss loci.
#'
#' @param loci a list of locus-specific sequence alignments of class `"matrix"` or `"DNAbin"`.
#' @param subset individual identifiers.
#' @param allele a regular expression defining an allele identifier (allows to distinguish individual IDs).
#' @param invert logical, whether to remove rather than retain the specified sequences.
#' @returns A numeric vector of locus-specific occupancies.
#' @export

subset_indiv <- function(loci, subset, allele=NULL, invert=FALSE) {

	is_matrix <- inherits(loci, "matrix")
	nloci <- ifelse(is_matrix, ncol(loci), length(loci))
	if (is.factor(subset)) {
		subset <- as.character(subset)
	}

	if (is_matrix) {
		part <- attr(loci, "part")
		if (is.character(subset)) {
			nams <- rownames(loci)
			if (isTRUE(invert)) {
				subset <- setdiff(nams, subset)
			}
			if (!is.null(allele)) {
				nams <- gsub(allele, "", nams)
				subset <- gsub(allele, "", subset)
			}
			loci <- loci[nams %in% subset,,drop=FALSE]	
		} else {
			loci <- loci[subset,,drop=FALSE]
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
			if (isTRUE(invert)) {
				nams <- rownames(loci[[i]])[!nams %in% subset]
			} else {
				nams <- rownames(loci[[i]])[nams %in% subset]
			}
			origattrib <- attributes(loci[[i]])
			loci[[i]] <- loci[[i]][nams,,drop=FALSE]
			allattrib <- c(attributes(loci[[i]]), origattrib)
			attributes(loci[[i]]) <- allattrib[match(names(origattrib), names(allattrib))]
		}
		totallattrib <- c(attributes(loci), totorigattrib)
		attributes(loci) <- totallattrib[match(names(totorigattrib), names(totallattrib))]
	}
	
	return(loci)	
}
