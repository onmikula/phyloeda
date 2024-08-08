#' Collapsed diploid genotypes.
#' 
#' @description
#' Collapses diploid genotypes into a single sequence using ambiguity codes.
#' 
#' @param loci an alignment in `matrix` or `DNAbin` format or a list of them.
#' @param allele regular expression specifying allele indicator.
#' @returns A modified alignment in `matrix` or `DNAbin` format or a list of them.
#' @export

collapse_diploid_genotypes <- function(loci, allele="_[012]$") {
	ambiguity_code <- c(AA="A", CC="C", GG="G", TT="T", "A-"="A", "C-"="C", "G-"="G", "T-"="T", "-A"="A", "-C"="C", "-G"="G", "-T"="T", AC="M", CA="M", AG="R", GA="R", AT="W", TA="W", CG="S", GC="S", CT="Y", TC="Y", GT="K", TG="K", "NN"="N", "N-"="N", "-N"="N", "--"="-")

	collapse <- function(x, allele, ambig) { 
		phased <- grepl(allele, rownames(x))
		nind <- sum(phased) / 2
		if (nind >= 1) {
			unphased <- x[!phased,,drop=FALSE]
			phased <- toupper(x[phased,,drop=FALSE])
			phased <- phased[order(rownames(phased)),,drop=FALSE]
			phased[!phased %in% c("A","C","G","T")] <- "-"
			indnames <- sub(allele, "", rownames(phased)[2*seq(nind)])
			phased <- setNames(lapply(seq(nind), function(i) unname(ambig[apply(phased[(2 * i) + c(-1, 0),,drop=FALSE], 2, paste, collapse="")])), indnames)
			phased <- do.call(rbind, phased)
			return(rbind(unphased, phased))
		} else {
			return(x)
		}
	}
	
	if (is.list(loci)) {
		loci <- 	lapply(loci, collapse, allele=allele, ambig=ambiguity_code)
	} else {
		loci <- collapse(loci, allele=allele, ambig=ambiguity_code)
	}
	return(loci)

}

