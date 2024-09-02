#' Writing sequence data.
#' 
#' @description
#' The functions writing sequence data into various file formats.
#'
#' `write_seq` is a wrapper for the following functions.
#'
#' `write_nexus` into .nexus file, with loci defined by an appropriate nexus block.
#'
#' `write_fasta` into .fasta file, loci concatenated or written into different files.
#'
#' `write_phylip` into .phy file, loci concatenated or written into different files.
#'
#' `write_ipyrad` into .alleles file used by ipyrad RADseq assembler.
#' 
#' @param loci a list of locus-specific sequence alignments of class `"matrix"` or `"DNAbin"`.
#' @param file character string, the name(s) of the file(s).
#' @param format character, the format of the file(s); it can be `"ipyrad"`, `"nexus"`, `"fasta"` or `"phylip"`.
#' @param nloci the number of loci to be written (the default `-1L` means all of them)
#' @param part either a two column matrix indicating starts and ends of loci in the concatenated aligments
#'   (in .fasta or .nexus files) or a name of partition file formatted according to RAxML or MrBayes convention.
#' @param rm.empty logical, whether to remove empty (e.g. all "-" or all "N") sequences.
#' @param interleaved logical, whether to use interleaved format; relevant for [write_nexus] and [write_phylip].
#' @param toupper logical, whether to represent nucleotides by upper case letters.
#' @param snps logical, whether to (extract and) include information about SNPs.
#' @param missing a symbol for missing data.
#' @param gap a symbol for gap in sequence.
#' @details a known issue is a correct specification of loci placed on negative strands in [write_ipyrad]
#'   applied to data from [read_stacks]. Currently, the functions do not distnguish missing data from gaps,
#'   although they use an appropriate coding convention (`missing="N"` and `gap="-"`). What is literally missing
#'   in the input is replaced by `"N"`, but dashes in input are left as they are, as the accepted inputs do not
#'   contain any information about their actual meaning.
#' @export

write_seq <- function(loci, file, format, nloci=-1L, part=NULL, rm.empty=FALSE, interleaved=FALSE, toupper=TRUE, snps=TRUE, missing="N", gap="-") {
	fun <- match.fun(paste("write", format, sep="_"))
	fun(loci=loci, file=file, nloci=nloci, part=part, rm.empty=rm.empty, interleaved=interleaved, toupper=toupper, snps=snps, missing=missing, gap=gap)
}



#' Write into ipyrad file.
#' 
#' @description
#' Writing sequence data into ipyrad .alleles file.
#' 
#' @param loci a list of locus-specific sequence alignments of class `"matrix"` or `"DNAbin"`.
#' @param file character string, the name of the .alleles file.
#' @param nloci the number of loci to be written (the default `-1L` means all of them)
#' @param part either a two column matrix indicating starts and ends of loci in the concatenated aligments
#'  (in .fasta or .nexus files) or a name of partition file formatted according to RAxML or MrBayes convention.
#' @param rm.empty logical, whether to remove empty (e.g. all "-" or all "N") sequences.
#' @param interleaved not relevant, just for compatibility with [read_seq].
#' @param toupper logical, whether to represent nucleotides by upper case letters.
#' @param snps logical, whether to (extract and) include information about SNPs.
#' @param missing a symbol for missing data.
#' @param gap a symbol for gap in sequence.
#' @export

write_ipyrad <- function(loci, file, nloci=-1L, part=NULL, rm.empty=FALSE, interleaved, toupper=TRUE, snps=TRUE, missing="N", gap="-") {

	if (!is.list(loci) | inherits(loci, "DNAbin")) {
		loci <- list(sequence=loci)
	}	
	if (isTRUE(rm.empty)) {
		loci <- remove_empty(loci, ambig=2)
	}

	if (isTRUE(snps)) {
		if (!is.null(attr(loci[[1]], "nalleles"))) {
			snps <- setNames(rep(list(snp=integer(0), pis=logical(0), nalleles=integer(0)), length(loci)), names(loci))
			for (i in seq_along(loci)) {
				snps[[i]]$snp <- seq(ncol(snps[[i]]))
				snps[[i]]$pis <- attr(snps[[i]], "pis")
				snps[[i]]$nalleles <- attr(snps[[i]], "nalleles")
			}	
		} else if (exists("find_snps")) {
			snps <- find_snps(loci)
		} else {
			snps <- FALSE
		}
	} else if (is.list(snps)) {
		if (!identical(names(snps[[1]]), c("snp","pis","nalleles"))) {
			snps <- FALSE
		}
	} else {
		snps <- FALSE
	}
	snpsstring <- lapply(sapply(loci, ncol), function(x) rep(" ", x))
	if (!isFALSE(snps)) {
		for (i in seq_along(loci)) {
			snpsstring[[i]][snps[[i]]$snp] <- ifelse(snps[[i]]$pis, "*", "-")
		}
	}
	snpsstring <- sapply(snpsstring, paste, collapse="")

	locino <- as.integer(sub("^[[:alpha:]]*", "", names(loci)))
	bedtable <- attr(loci, "bedtable")
	if (!is.null(bedtable)) {
		gencoor <- paste0(bedtable[,1], ":", bedtable[,2]+1, "-", bedtable[,3]+1)  # 'bedtable[,3]+1' only for compatibility with original ipyrad output
		snpsstring <- paste0(snpsstring, "|", locino, ":", gencoor, "|")
	} else {
		gencoor <- rep("", length(loci))
		snpsstring <- paste0(snpsstring, "|", locino, "|")
	}
	
	if (inherits(loci[[1]], "DNAbin")) {
		loci <- lapply(loci, as.matrix)
	}
	if (isTRUE(toupper)) {
		loci <- lapply(loci, toupper)
	}
	seqnam <- sort(unique(unlist(lapply(loci, rownames))))
	nchnam <- c(5, setNames(nchar(seqnam), seqnam))
	offset <- ceiling(0.2 * max(nchnam + 5)) / 0.2
	seqoff <- sapply(offset - nchnam[-1], function(x) paste(rep(" ", x), collapse=""))
	for (i in seq_along(loci)) {
		nams <- mapply(paste0, rownames(loci[[i]]), seqoff[rownames(loci[[i]])], SIMPLIFY=TRUE)
		loci[[i]] <- paste0(nams, apply(loci[[i]], 1, paste, collapse=""))
	}
	snpsstring <- paste0("//", paste0(rep(" ", offset - 2)), snpsstring)
	if (!isFALSE(snps)) {
		loci <- Map(c, loci, snpsstring)
	}
	
	loci <- unlist(loci)
	writeLines(loci, con=file, sep="\n", useBytes=FALSE)
	
}



#' Write into nexus file.
#' 
#' @description
#' Writing sequence data into nexus file.
#' 
#' @param loci a list of locus-specific sequence alignments of class `"matrix"` or `"DNAbin"`.
#' @param file character string, the name of the .alleles file.
#' @param nloci the number of loci to be written (the default `-1L` means all of them)
#' @param part either a two column matrix indicating starts and ends of loci in the concatenated aligments
#'  (in .fasta or .nexus files) or a name of partition file formatted according to RAxML or MrBayes convention.
#' @param rm.empty logical, whether to remove empty (e.g. all "-" or all "N") sequences.
#' @param interleaved logical, whether to use interleaved format.
#' @param toupper logical, whether to represent nucleotides by upper case letters.
#' @param snps logical, whether to (extract and) include information about SNPs.
#' @param missing a symbol for missing data.
#' @param gap a symbol for gap in sequence.
#' @export

write_nexus <- function(loci, file, nloci=-1L, part=NULL, rm.empty=FALSE, interleaved=FALSE, toupper=TRUE, snps=TRUE, missing="N", gap="-") {

	bedtable <- attr(loci, "bedtable")
	if (is.null(part)) {
		part <- attr(loci, "part")
	}
	if (!is.list(loci) | inherits(loci, "DNAbin")) {
		loci <- list(sequence=loci)
	}
	if (isTRUE(rm.empty)) {
		loci <- remove_empty(loci, ambig=2)
	}
	if (!is.null(attr(loci[[1]], "snp"))) {
		snpattr <- lapply(loci, function(x) cbind(snp=attr(x, "snp"), pis=as.character(attr(x, "pis")), nalleles=attr(x, "nalleles")))
		snpattr <- cbind(locus=rep(names(loci), sapply(snpattr, nrow)), do.call(rbind, snpattr))
	} else {
		snpattr <- NULL
	}

	nbp <- sapply(loci, ncol)
	nz <- nbp > 0
	nbp <- nbp[nz]
	loci <- loci[nz]
	seqnames <- sort(unique(unlist(lapply(loci, rownames))))
	ntaxa <- length(seqnames)
	nsites <- sum(nbp)
	if (!is.null(part)) {
		part <- part[nz,,drop=FALSE]
	} else {
		cs <- cumsum(nbp)
		part <- matrix(c(1, cs[-length(nbp)] + 1, cs), length(nbp), 2, dimnames=list(names(nbp), NULL))
	}
	if (inherits(loci[[1]], "matrix")) {
		empty <- lapply(nbp, function(x) rep(missing, x))
	} else if (inherits(loci[[1]], "DNAbin")) {
		empty <- lapply(nbp, function(x) ape::as.DNAbin(rep(missing, x)))
	}

	if (inherits(loci[[1]], "DNAbin")) {
		loci <- lapply(loci, as.matrix)
	}
	if (isTRUE(toupper)) {
		loci <- lapply(loci, toupper)
	}
	sequences <- setNames(vector("list", ntaxa), seqnames)
	for (i in seq(ntaxa)) {
		matches <- sapply(lapply(loci, rownames), function(x) match(seqnames[i], x))
		sequences[[i]] <- lapply(seq_along(loci), function(j) if (!is.na(matches[j])) loci[[j]][matches[j],] else character(0))
		loci <- lapply(seq_along(loci), function(j) if (!is.na(matches[j])) loci[[j]][-matches[j],,drop=FALSE] else loci[[j]])
		zero <- sapply(sequences[[i]], length) == 0
		sequences[[i]][zero] <- empty[zero]
		sequences[[i]] <- unlist(sequences[[i]])
	}
	if (isTRUE(interleaved)) {
		padded <- nchar(seqnames)
		padded <- max(padded) - padded + 4
		padded <- sapply(lapply(padded, function(times) rep(" ", times)), paste, collapse="")
		padded <- paste0(seqnames, padded)
		for (i in seq(ntaxa)) {
			len <- length(sequences[[i]])
			lines <- findInterval(1:len, vec=seq(1, len, by=80), rightmost.closed=FALSE, all.inside=FALSE, left.open=FALSE)
			sequences[[i]] <- unname(sapply(Map(c, padded[i], split(sequences[[i]], lines)), paste, collapse=""))
		}
		nlines <- length(sequences[[i]])
		sequences <- unname(unlist(sequences)[rep(seq(nlines), each=ntaxa) + rep(nlines * (seq(ntaxa)-1), nlines)])
		blocks <- which(seq_along(sequences) %% ntaxa == 0)[-nlines]
		sequences[blocks] <- paste0(sequences[blocks], "\n")
	} else {
		sequences <- lapply(sequences, paste, collapse="")
		sequences <- unname(c(seqnames, unlist(sequences))[rep(seq(ntaxa), each=2) + (seq(2 * ntaxa) %% 2 == 0) * ntaxa])
	}

	binary <- any(loci[[1]] %in% c(0, 1, 2))
	if (isFALSE(binary)) {
		format <- paste("Format", "datatype=DNA", paste0("gap=", gap), paste0("missing=", missing, ";"))
	} else {
		format <- paste("Format", "datatype=integerdata", symbols="\"012\"", paste0("gap=", gap), paste0("missing =", missing, ";"))
	}
	header <- c("#NEXUS\n", "Begin data;", paste0("dimensions", "  ", "ntax=", ntaxa, " ", "nchar=", nsites, ";"), format, "Matrix\n")

	nexus <- c(header, sequences, ";\n", "End;")
	if (!is.null(part)) {
		p <- apply(part, 1, function(x) paste(x, collapse=" - "))
		p <- paste(names(p), p, sep=" = ")
		p <- paste(paste(paste("\t", "charset", sep=""), p), ";", sep="")
		p <- c("Begin assumptions;", p, "End;")
		nexus <- c(nexus, "\n", p)
	}

	if (!is.null(bedtable)) {
		bedtable <- bedtable[nz,,drop=FALSE]
		bed <- c("Begin bedtable;", gsub(" ", "", apply(bedtable, 1, paste, collapse="\t")), "End;")
		nexus <- c(nexus, "\n", bed)
	}	
	if (!is.null(snpattr)) {
		snpattr <- c("Begin snps;", paste(colnames(snpattr), collapse="\t"), gsub(" ", "", apply(snpattr, 1, paste, collapse="\t")), "End;")
		nexus <- c(nexus, "\n", snpattr)
	}

	writeLines(nexus, con=file, sep="\n", useBytes=FALSE)

}


