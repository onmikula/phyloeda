#' Sequence matches.
#' 
#' @description
#' Classify sequences into taxa based on their identity to already classified ones.
#'
#' @param taxa a data frame with information about classification of haplotypes to taxa (species, clades, lineages)
#'   or a name of file with such classification. Haplotype names are in the 1st column, taxa in the 2nd column.
#' @param haps either an output of [haplotypes] or a data frame with assignment of sequences to haplotypes
#'   (`assign` component of the list produced by [haplotypes] or a name of tab-delimited file with such data).
#' @param haps either a data frame with assignment of sequences to haplotypes (as produced by [haplotypes])
#'   or a name of file with such an assignment. Sequence names are in the 1st column, haplotype names in the 2nd column.
#' @param na.rm logical, whether to omit sequences that could not be classified into taxa.
#' @param retain logical, whether to retain all columns of `taxa` data frame or just the first two of them.
#' @param col.names optional character vector, column names of the output data frame.
#' @returns A data frame with similar to `taxa` input, but having in the 1st column sequence names
#'   and in the 2nd column assignment to taxa, based on the correspondence of haplotype names.
#' @export

seqmatches <- function(taxa, haps, na.rm=FALSE, retain=TRUE, col.names) {
	if (is.character(taxa) & length(taxa) == 1) {
		taxa <- read.delim(taxa)[,1:2]
	} else {
		taxa <- as.data.frame(taxa)
	}
	if (is.character(haps) & length(haps) == 1) {
		haps <- read.delim(haps, col.names=c("seq", "hap"))[,1:2]
	} else if (is.list(haps) & !is.data.frame(haps)) {
		haps <- haps$assign
	}
	haps <- haps[!haps$seq %in% taxa[,1],,drop=FALSE]
	if (nrow(haps) > 0) {
		matching <- match(haps$hap, taxa[,1])
		result <- data.frame(seq=haps$seq, tax=taxa[matching,2])	
		if (isTRUE(retain) & length(taxa) > 2) {
			result <- data.frame(result, taxa[matching,3:length(taxa)])
		}
	} else {
		result <- data.frame(seq=character(0), tax=character(0))
		if (isTRUE(retain) & length(taxa) > 2) {
			result <- data.frame(result, taxa[FALSE,3:length(taxa)])
		}
	}
	if (!missing(col.names)) {
		names(result)[1:length(col.names)] <- col.names
	}
	if (isTRUE(na.rm)) {
		result <- result[!is.na(result[,2]),]
	}
	return(result)
}


#' Query sequences.
#' 
#' @description
#' Finds query sequences for phylogenetic placement, i.e., those that are different from reference sequences and from each other.
#'
#' @param ref a character vector with IDs of reference sequences or a single tree (readable by [read_tree]),
#'   whose tip labels are used as IDs of reference sequences.
#' @param haps either an output of [haplotypes] or a data frame with assignment of sequences to haplotypes
#'   (`assign` component of the list produced by [haplotypes] or a name of tab-delimited file with such data).
#' @returns A character vector with IDs of query sequences.
#' @export

findqueryseq <- function(ref, haps) {
	if ((is.character(ref) & length(ref) == 1) | inherits(ref, "phylo")) {
		ref <- read_tree(ref)$tip.label
	} 
	if (is.character(haps) & length(haps) == 1) {
		haps <- read.delim(haps, col.names=c("hap", "seq"))[,1:2]
	} else if (is.list(haps) & !is.data.frame(haps)) {
		haps <- haps$assign
	}
	haps$query <- TRUE
	haps <- split(haps, haps$hap)
	for (i in seq_along(haps)) {
		if (any(haps[[i]]$seq %in% ref)) {
			haps[[i]]$query <- FALSE
		} else {
			haps[[i]]$query <- rep(c(TRUE, FALSE), c(1, nrow(haps[[i]]) - 1))
		}
	}
	haps <- do.call(rbind, haps)
	query <- sort(haps$seq[haps$query])
	return(query)
}


#' EPA input.
#' 
#' @description
#' Prepares sequence file to be used as an input for phylogenetic placement, i.e., including reference and query sequences.
#'
#' @param ref a character vector with IDs of reference sequences or a single tree (readable by [read_tree]),
#'   whose tip labels are used as IDs of reference sequences.
#' @param query a character vector with IDs of query sequences. Also it could be an object accepted as `haps` argument
#'   by [findqueryseq], i.e., either an output of [haplotypes] or a data frame with assignment of sequences to haplotypes
#'   (`assign` component of the list produced by [haplotypes] or a name of tab-delimited file with such data).
#' @param seq either a a sequence alignment or a name of .fasta file which contains both reference and query sequences.
#' @param file character string, the name of .fasta file with the resulting alignment.
#' @param return logical, whether to return the resulting sequence alignment.
#' @returns If `return=TRUE`, the resulting sequence alignment in the same format in which it was supplied or imported in.
#' @export

epainput <- function(ref, query, seq, file, return=FALSE) {
	if ((is.character(seq) & length(seq) == 1)) {
		seq <- read.fasta(seq)
	}
	if ((is.character(ref) & length(ref) == 1) | inherits(ref, "phylo")) {
		ref <- read_tree(ref)$tip.label
	} 
	if (any(!ref %in% rownames(seq))) {
		stop("At least one of the reference sequences is not in the alignment.")
	}
	if (is.character(query) & length(query) == 1) {
		if (file.exists(query)) {
			query <- read.delim(query, col.names=c("hap", "seq"))[,1:2]
		}
	} else if (is.list(query) & !is.data.frame(query)) {
		query <- query$assign
	}
	if (is.data.frame(query) | is.matrix(query)) {
		query <- findqueryseq(ref, haps=query)
	}
	query <- intersect(query, rownames(seq))
	if (length(query) == 0) {
		warning("None of the query sequences is in the alignment, NULL is returned.")
		return(NULL)
	} else {
		epa <- seq[c(ref, query),]
		if (!missing(file)) {
			write.fasta(epa, file)
		}
		if (isTRUE(return)) {
			return(epa)
		}	
	}
}


#' EPA format.
#' 
#' @description
#' Summarizes classification of sequences into taxa in a format similar to output of [epaclades::best_placement()].
#'
#' @param taxa a data frame with information about classification of haplotypes to taxa (species, clades, lineages)
#'   or a name of file with such classification. Haplotype names are in the 1st column, taxa in the 2nd column.
#' @param haps either an output of [haplotypes] or a data frame with assignment of sequences to haplotypes
#'   (`assign` component of the list produced by [haplotypes] or a name of tab-delimited file with such data).
#'   Sequence names are in the 1st column, haplotype names in the 2nd column.
#' @param epa an output of [epaclades::classify_jplace()], [epaclades::classify_sequences()] or
#'   [epaclades::best_placement()] or a name of tab-delimited file with output of [epaclades::best_placement()].
#'   If missing, the epa formatting is still imposed to the output based on `taxa` and `haps`.
#' @param file character string, the name of the output file.
#' @param return logical, whether to return the resulting sequence alignment.
#' @param na.rm logical, whether to omit sequences that could not be classified into taxa.
#' @param method optional character string, method of classification into units predefined in `taxa`.
#' @returns A data frame with six columns labelled "id" (sequence name), "clade" (taxon), "method",
#'   "position" (stem or crown), "probability" (for the position), "totalprob" (summed over stem and crown branches).
#'   The last three columns are relevant only for classification by phylogenetic placement.
#'   The "method" column takes three possible values: "predefined" (or the name of the original delimitation method
#'   from `method` argument), "epa" (classification by phylogenetic placement) or "match" (if the classification
#'   was inherited from another through sequence matching, possibly with the associated probabilities).
#' @export

epaformat <- function(taxa, haps, epa, file, return=TRUE, na.rm=FALSE, method) {
	taxa <- as.data.frame(taxa)
	delim <- grep("^method", names(taxa), ignore.case=TRUE, value=TRUE)
	if (length(delim) == 1) {
		taxa <- setNames(data.frame(taxa[,1:2], taxa[,delim]), c("id", "clade", "method"))
	} else if (!missing(method)) {
		taxa <- setNames(data.frame(taxa[,1:2], method), c("id", "clade", "method"))
	} else {
		taxa <- setNames(data.frame(taxa[,1:2], "predefined"), c("id", "clade", "method"))
	}
	taxa$totalprob <- taxa$probability <- taxa$position <- NA		

	if (missing(epa)) {
		epa <- NULL
	}
	if (!is.null(epa)) {
		if (is.character(epa) & length(haps) == 1) {
			epa <- setNames(read.delim(epa)[,1:5], c("id", "clade", "position", "probability", "totalprob"))
		} else if (inherits(epa, "jplace")) {
			epa <- epaclades::best_placement(epaclades::classify_sequences(epa))
			epa <- epa[, c("id", "clade", "position", "probability", "totalprob")]
		} else if (is.data.frame(epa)) {
			if (!"totalprob" %in% names(epa)) {
				epa <- epaclades::best_placement(epa)
			}
		}
		epa$method <- "epa"
		epa <- epa[,c("id", "clade", "method", "position", "probability", "totalprob")]
		result <- rbind(taxa, epa)
	} else {
		result <- taxa
	}

	if (is.character(haps) & length(haps) == 1) {
		haps <- setNames(read.delim(haps)[,1:2], c("seq", "hap"))
	} else if (is.list(haps) & !is.data.frame(haps)) {
		haps <- haps$assign
	} else if (is.data.frame(haps)) {
		haps <- setNames(haps[,1:2], c("seq", "hap"))
	}
	haps <- seqmatches(taxa=result, haps=haps, na.rm=na.rm, retain=TRUE, col.names=c("id", "clade"))
	if (nrow(haps) > 0) {
		haps$method <- "match"
	}	
	haps <- haps[,c("id", "clade", "method", "position", "probability", "totalprob")]
	result <- rbind(result, haps)
	rownames(result) <- NULL

	if (!missing(file)) {
		write.delim(result, file)
	}
	if (isTRUE(return)) {
		return(result)
	}
	
}

