#' Unique haplotypes.
#' 
#' @description
#' Extracts unique haplotypes from an alignment.
#'
#' @param input a character string with file name (fasta or phylip interleaved format)
#'   or sequence alignment in `DNAbin` or `matrix` format.
#' @param fasta character, a name of the output .fasta file with the alignment of haplotypes.
#' @param txt character, a name of the output .txt file with assignment of sequences to haplotypes.
#' @param return logical, whether to return the output as a list.
#' @param outgroup a character vector with outgroup sequence names or their unique identifier(s),
#'   e.g., a character string `"outgroup"`. The outgroups are excluded from haplotype identification,
#'   but presented in the output.
#' @param minpres numeric, a minimum required proportion of nucleotides non-missing in the sequence 
#'   under consideration and also in a haplotype it could be assigned to (default is 1).
#' @param maxmiss numeric, a maximum acceptable proportion of nucleotides missing in the sequence
#'   under consideration, but non-missing in a haplotype it could be assigned to (default is 1).
#' @details Under the default values of `minpres` and `maxmiss` a shorter sequence is considered identical
#'   to a longer one, if it does not contain any nucleotides missing in the longer sequence, but not vice versa.
#'   If `minpres=1` and `maxmiss=0`, a precise match is required.
#' @returns If `return == TRUE`, a list with components `haplotypes` (alignment of haplotypes),
#'   `assign` (assignment of sequences to haplotypes), `no` (number of unique haplotypes found)
#'   and `outgroup` (alignment of outgroup sequences).
#' @export

haplotypes <- function(input, fasta=NULL, txt=NULL, return=TRUE, outgroup=NULL, minpres=1, maxmiss=1) {
	if (inherits(input, "DNAbin")) {
		seq <- input <- toupper(as.character(input))
	} else if (inherits(input, "matrix")) {
		seq <- input <- toupper(input)
	} else {
		if (substr(readLines(input, n=1), 1, 1) == ">") {
			seq <- input <- toupper(ape::read.dna(input, format="fasta", as.character=TRUE, as.matrix=TRUE))
		} else {
			seq <- input <- toupper(ape::read.dna(input, as.character=TRUE, as.matrix=TRUE))
		}	
	}
	if (!is.null(outgroup)) {
		o <- unlist(lapply(outgroup, grep, x=rownames(seq)))
		o <- sort(unique(c(o, unlist(lapply(outgroup, grep, x=rownames(seq), fixed=TRUE)))))
		outgr <- seq[o,,drop=FALSE]
		seq <- seq[-o,,drop=FALSE]
		if (length(outgr) == 0) {
			outgroup <- NULL
			warning("alignment doesn't contain any of the specified outgroups")
		} else {
			outgroup <- outgr
		}
	}
	seq <- matrix(ifelse(seq %in% c("A","C","T","G"), seq, NA), nrow(seq), ncol(seq), dimnames=dimnames(seq))
	seq <- seq[order(rowSums(!is.na(seq)), decreasing=TRUE),]
	dat <- !is.na(seq)
	cvr <- rep(1, nrow(seq)) %*% t(rowSums(dat))
	ovr <- dat %*% t(dat)
	del <- (t(cvr) - ovr) / t(cvr)
	ovr <- ovr / cvr
	seq[is.na(seq)] <- "-"
	dst <- ape::dist.dna(x=ape::as.DNAbin(seq), model="raw", pairwise.deletion=TRUE, as.matrix=TRUE)
	if (all(dst[upper.tri(dst)] > 0)) {
		haps <- setNames(vector("list", nrow(dst)), rownames(dst))
	} else {
		haps <- as.list(apply(dst == 0 & ovr >= minpres & del <= maxmiss & upper.tri(dst), 1, which))
		haps <- lapply(haps, names)
		for (i in seq_along(haps)[-1]) {
			haps[[i]] <- setdiff(haps[[i]], unlist(haps[1:(i-1)]))
		}
		haps <- haps[!names(haps) %in% unlist(haps)]
	}
	haps <- Map(c, names(haps), haps)
	haps <- data.frame(seq=unlist(haps), hap=rep(names(haps), sapply(haps, length)), row.names=NULL)
	seq <- input[rownames(input) %in% haps$hap,,drop=FALSE]
	seq[is.na(seq)] <- "-"
	if (!is.null(outgroup)) {
		seq <- rbind(seq, outgroup)
	}
#	if (!is.null(fasta) & is.null(txt)) txt <- sub("\\.[[:alpha:]]+$", ".txt", fasta)
#	if (!is.null(txt) & is.null(fasta)) fasta <- sub("\\.[[:alpha:]]+$", ".fasta", txt)
	if (!is.null(fasta)) {
		ape::write.dna(x=seq, file=fasta, format="fasta")
	}
	if (!is.null(txt)) {
		utils::write.table(haps, file=txt, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
	}
	if (isTRUE(return)) {
		return(list(haplotypes=seq, assign=haps, no=nrow(seq), outgroup=outgroup))
	}
}



#' Query sequences.
#' 
#' @description
#' Finds query sequences for phylogenetic placement, i.e., those that are different from reference sequences and from each other.
#'
#' @param ref a character vector with IDs of reference sequences or a single tree (readable by [read_tree]),
#'   whose tip labels are used as IDs of reference sequences.
#' @param haps either an output of [haplotypes] or a data frame with assignment of sequences to haplotypes (`assign` component
#'   of the list produced by [haplotypes] or a name of tab-delimited file with such data.
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



#' Classify sequences.
#' 
#' @description
#' Classify sequences based on their identity to already classified ones.
#'
#' @param info a data frame with information about classification of haplotypes to taxa (species, clades, lineages)
#'   or a name of file with such classification. Haplotype names are in the 1st column, taxa in the 2nd column.
#' @param haps either a data frame with assignment of sequences to haplotypes (as produced by [haplotypes])
#'   or a name of file with such an assignment. Haplotype names are in the 1st column, sequence names in the 2nd column.
#' @param col.names column names of the output data frame.
#' @param na.rm logical, whether to omit sequences whose classification to taxa cannot be established
#'   from the classification of haplotypes.
#' @returns A data frame with similar to `info` input, but having in the 1st column sequence names
#'   and in the 2nd column assignment to taxa, based on the correspondence of haplotype names.
#' @export

classifyseq <- function(info, haps, col.names=c("seq", "tax"), na.rm=FALSE) {
	if (is.character(info) & length(info) == 1) {
		info <- read.delim(info, col.names=c("hap", "tax"))[,1:2]
	}
	if (is.character(haps) & length(haps) == 1) {
		haps <- read.delim(haps, col.names=c("hap", "seq"))[,1:2]
	}
	result <- data.frame(haps$seq, info$tax[match(haps$hap, info$hap)], col.names=col.names)
	if (isTRUE(na.rm)) {
		result <- result[!is.na(result[,2]),]
	}
	return(result)
}


#' EPA input.
#' 
#' @description
#' Prepares sequence file to be used as an input for phylogenetic placement, i.e., including reference and query sequences.
#'
#' @param ref a character vector with IDs of reference sequences or a single tree (readable by [read_tree]),
#'   whose tip labels are used as IDs of reference sequences.
#' @param query a character vector with IDs of query sequences. Also it could be an object accepted as `haps` argument by [findqueryseq],
#'   i.e., either an output of [haplotypes] or a data frame with assignment of sequences to haplotypes (`assign` component
#'   of the list produced by [haplotypes] or a name of tab-delimited file with such data.
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
		stop("None of the query sequences is in the alignment.")
	}
	epa <- seq[c(ref, query),]
	if (!missing(file)) {
		write.fasta(epa, file)
	}
	if (isTRUE(return)) {
		return(epa)
	}
}



