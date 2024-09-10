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
		if (length(o) == 0) {
			outgroup <- NULL
			warning("alignment doesn't contain any of the specified outgroups")
		} else {
			outgroup <- seq[o,,drop=FALSE]
			seq <- seq[-o,,drop=FALSE]
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
