#' Most divergent sequences.
#' 
#' @description
#' Selection of a specified number of sequences, which are most divergent from each other.
#'
#' @param input a character string with name of the file in format accepted by [read_seq]
#'   or a sequence alignment in `DNAbin` or `matrix` format or an output of [haplotypes].
#' @param n number of sequences to retain.
#' @param minlen numeric the minimum required length of sequence.
#' @param fasta character, a name of the output .fasta file with the alignment of haplotypes.
#' @param return logical, whether to return the output as a list.
#' @param outgroup a character vector with outgroup sequence names or their unique identifier(s),
#'   e.g., a character string `"outgroup"`.
#' @param reserved a character vector listing sequences that must be included in the selection.
#' @param model a character string specifying the nucleotide substitution model used to correct for multiple
#'   mutations when calculating pairwise distences between sequences; deafult is `"raw"` (no correction).
#' @details The algorithm proceeds by stepwise elimination of one of the most similar sequences in the alignment.
#' @returns An alignment of `n` retained sequences (in `DNAbin` or `matrix` format, depending on the input).
#' @export

divergentseq <- function(input, n, minlen=0, fasta=NULL, return=TRUE, outgroup=NULL, reserved=NULL, model="raw") {
	if (inherits(input, "DNAbin")) {
		seq <- toupper(as.character(input))
	} else if (inherits(input, "matrix")) {
		seq <- toupper(input)
	} else if (is.list(input) & isTRUE(all(c("haplotypes", "assign") %in% names(input)))) {
		seq <- input <- input$haplotypes
	} else {
		if (substr(readLines(input, n=1), 1, 1) == ">") {
			input <- ape::read.dna(input, format="fasta", as.character=TRUE, as.matrix=TRUE)
			seq <- toupper(input)
		} else {
			input <- ape::read.dna(input, as.character=TRUE, as.matrix=TRUE)
			seq <- toupper(input)
		}	
	}
	if (!is.null(outgroup)) {
		o <- unlist(lapply(outgroup, grep, x=rownames(seq)))
		o <- sort(unique(c(o, unlist(lapply(outgroup, grep, x=rownames(seq), fixed=TRUE)))))
		outgr <- rownames(seq)[o]
		seq <- seq[-o,,drop=FALSE]
	} else {
		outgr <- NULL
	}
	if (any(!reserved %in% rownames(seq))) {
		reserved <- intersect(reserved, rownames(seq))
		if (length(reserved) == 0) {
			reserved <- NULL
			warning("the 'reserved' sequences are not present in the alignment")
		} else {
			warning("some of the 'reserved' sequences are not present in the alignment and were removed from the list")
		}
	}
	if (length(reserved) > n) {
		stop("the reserved list is longer than the number of sequences to be selected")
	}
	res <- rownames(seq) %in% reserved	
	len <- seqlen(seq)
	ord <- order(res, len, decreasing=TRUE)
	seq <- seq[ord,,drop=FALSE]
	seq <- seq[len[ord] >= minlen | res[ord],,drop=FALSE]
	dst <- ape::dist.dna(ape::as.DNAbin(seq), pairwise.deletion=TRUE, model=model, as.matrix=TRUE)
	dst[reserved, reserved] <- NA
	while (nrow(dst) > n) {
		mindst <- min(dst[upper.tri(dst)], na.rm=TRUE)
		ij <- which(dst == mindst & upper.tri(dst), arr.ind=TRUE)
		excl <- ij[order(ij[,1], ij[,2]),,drop=FALSE][1,2]
		dst <- dst[-excl, -excl]
	}
	seq <- input[rownames(input) %in% c(rownames(dst), outgr),,drop=FALSE]
	if (!is.null(fasta)) {
		ape::write.dna(x=seq, file=fasta, format="fasta")
	}
	if (isTRUE(return)) {
		return(seq)
	}
}



#' Most distant cases.
#' 
#' @description
#' Selection of a specified number of cases (e.g., sequences or tree tips), which are most distant from each other.
#'
#' @param dst distance matrix, object of class `dist` or `matrix`, the latter is assumed to be square and symmetric.
#' @param n number of cases to retain.
#' @param phy an object of class `phylo` with phylogeny to be subset.
#' @param seq an object of class `DNAbin` or `matrix` with sequence alignment to be subset.
#' @param outgroup a character vector with outgroup names or their unique identifier(s),
#' @param reserved A character vector listing sequences that must be included in the selection.
#' @param model A character string specifying the nucleotide substitution model used to correct for multiple
#'   mutations when calculating pairwise distences between sequences; deafult is `"raw"` (no correction).
#' @details The algorithm operates on the supplied distance matrix. It can be used to subsample tree tips based on distances
#'   between them (if `dst=cophenetic(phy)`) or sequences based on their mutual distance (if `dst=dist.dna(seq)`).
#'   However, phylogenetic distances can be also used to subsample sequences and vice versa, or the distance matrice can be
#'   of entirely different kind, e.g., matrix of geographical distances between individuals or genetic distances calculated
#'   from another part of genome. Anyway, either `phy` or `seq` or none of them, but not both at the same time must be supplied.
#' @returns If `phy` or `seq` are supplied the function return subset phylogenetic tree or sequence alignment, respectively.
#'   Otherwise, it returns a vector with IDs of cases to be retained.
#' @export

maxdist <- function(dst, n, phy=NULL, seq=NULL, outgroup=NULL, reserved=NULL, ord=NULL) {
	if (!is.null(phy) & !is.null(seq)) {
		stop("either 'phy' or 'seq' or none of them, but not both at the same time must be supplied")
	}
	dst <- as.matrix(dst)
	rownam <- rownames(dst)
	if (!is.null(outgroup)) {
		o <- unlist(lapply(outgroup, grep, x=rownames(dst)))
		o <- sort(unique(c(o, unlist(lapply(outgroup, grep, x=rownames(dst), fixed=TRUE)))))
		outgr <- rownames(dst)[o]
		dst <- dst[-o,-o,drop=FALSE]
	} else {
		outgr <- NULL
	}
	if (any(!reserved %in% rownames(dst))) {
		dst <- intersect(dst, rownames(dst))
		if (length(dst) == 0) {
			dst <- NULL
			warning("the 'reserved' IDs are not present in the distance matrix")
		} else {
			warning("some of the 'reserved' IDs are not present in the distance matrix and were removed from the list")
		}
	}
	if (length(reserved) > n) {
		stop("the reserved list is longer than the number of sequences to be selected")
	}
	res <- rownames(dst) %in% reserved
	if (is.null(ord)) {
		ord <- seq(nrow(dst))
	} else {
		ord <- ord[intersect(rownames(dst), names(ord))]
	}
	ord <- order(res, ord, decreasing=c(TRUE, FALSE), method="radix")
	dst <- dst[ord,ord]
	dst[reserved, reserved] <- NA
	while (nrow(dst) > n) {
		mindst <- min(dst[upper.tri(dst)], na.rm=TRUE)
		ij <- which(dst == mindst & upper.tri(dst), arr.ind=TRUE)
		excl <- ij[order(ij[,1], ij[,2]),,drop=FALSE][1,2]
		dst <- dst[-excl, -excl]
	}
	retain <- intersect(rownam, c(rownames(dst), outgr))
	if (!is.null(phy)) {
		rooted <- ape::is.rooted(phy)
		phy <- ape::drop.tip(phy, tip=setdiff(phy$tip.label, retain))
		if (isFALSE(rooted)) {
			phy <- ape::unroot(phy)
		}
		return(phy)
	} else if (!is.null(seq)) {
		seq <- seq[intersect(rownames(seq), retain),,drop=FALSE]
		return(seq)
	} else {
		return(retain)
	}
}

