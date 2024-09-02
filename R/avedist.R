#' Average genetic distances.
#' 
#' @description
#' Calculates average genetic distance between the specified species (populations, lineages) and within them.
#' 
#' @param loci The list of locus-specific alignments or a single such alignment (possibly 
#'   with concatenated multi-locus data). The row names must be interpretable as individual labels
#'   (possibly with help of `seqmap` or `allele` arguments).
#' @param info A matrix or data frame (1st column individuals, 2nd column species)
#'   or a name of a tab-delimited file with such data.
#' @param seqmap A data frame mapping names of sequences (1st column) to those of individuals (2nd column).
#' @param allele An allele identifier, regular expression distinguishing sequences from the same individual.
#' @param model The nucleotide substitution model indicated as in [ape::dist.dna()].
#' @param diag Logical, whether to calculate average intraspecific distances, default is `TRUE`.
#' @returns A square symmetric matrix giving mean genetic distance between species (off-diagonal)
#'   and within them (diagonal).
#' @export

avedist <- function(loci, info=NULL, seqmap=NULL, allele=NULL, model="raw", diag=TRUE) {

	if (is.list(loci)) {
		loci <- concatenate(loci, as.DNAbin=TRUE)
	} else if (!inherits(loci, "DNAbin")) {
		loci <- ape::as.DNAbin(loci)
	}

	im <- !is.null(info)
	sm <- !is.null(seqmap)
	al <- !is.null(allele)
	if (im) {
		if (is.character(info) & length(info) == 1) {
			info <- read.delim(info, col.names=c("Individual", "Species"))
		} else {
			info <- setNames(as.data.frame(info[,1:2]), c("Individual", "Species"))
		}
		if (!sm) {
			if (al) {
				seqmap <- data.frame(Sequence=paste0(info[,1], allele), Individual=info[,1])
			} else {
				seqmap <- data.frame(Sequence=info[,1], Individual=info[,1])
			}
		} else if (al) {
			seqmap <- setNames(as.data.frame(seqmap[,1:2]), c("Sequence", "Individual"))
			seqmap <- unique(rbind(seqmap, data.frame(Sequence=paste0(info[,1], allele), Individual=info[,1])))
			seqnames <- lapply(seqmap$Sequence, grep, x=rownames(loci), value=TRUE)
			seqmap <- unique(data.frame(Sequence=unlist(seqnames), Individual=rep(seqmap$Individual, sapply(seqnames, length))))
		}
	} else {
		if (!sm) {
			if (al) {
				indnames <- sub(allele, "", rownames(loci))
				info <- data.frame(Individual=unique(indnames), Species=unique(indnames))
				seqmap <- data.frame(Sequence=rownames(loci), Individual=indnames)
			} else {
				info <- data.frame(Individual=rownames(loci), Species=rownames(loci))
				seqmap <- data.frame(Sequence=rownames(loci), Individual=rownames(loci))
			}
		} else {
			seqmap <- setNames(as.data.frame(seqmap[,1:2]), c("Sequence", "Individual"))
			info <- unique(data.frame(Individual=seqmap$Individual, Species=seqmap$Individual))
		}
	}

	rows <- lapply(seqmap$Sequence, grep, x=rownames(loci))
	indnames <- seqmap$Individual[rep(seq_along(seqmap$Sequence), sapply(rows, length))][order(unlist(rows))]
	sp <- factor(info$Species[match(indnames, info$Individual)])
	dst <- ape::dist.dna(loci, model=model, pairwise.deletion=TRUE, as.matrix=TRUE)
	ave <- matrix(, nlevels(sp), nlevels(sp), dimnames=list(levels(sp), levels(sp)))
	for (i in seq(nlevels(sp) - 1)) {
		ii <- sp == levels(sp)[i]
		for (j in (i + 1):nlevels(sp)) {
			jj <- sp == levels(sp)[j]
			ave[i,j] <- ave[j,i] <- mean(dst[ii,jj], na.rm=TRUE)
		}
	}
	
	if (isTRUE(diag)) {
		for (i in seq(nlevels(sp))) {
			ii <- sp == levels(sp)[i]
			if (sum(ii) > 1) {
				ave[i,i] <- mean(dst[ii,ii][lower.tri(diag(sum(ii)))], na.rm=TRUE)
			}
		}
	}

	return(ave)

}
