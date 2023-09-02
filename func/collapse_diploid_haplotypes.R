collapse_diploid_haplotypes <- function(loci, allele="_[[:digit:]]$") {
#	iupac_ambiguity <- c(A="A", C="C", G="G", T="T", M="AC", R="AG", W="AT", S="CG", Y="CT", K="GT", V="ACG", H="ACT", D="AGT", B="CGT", N="ACGT")
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
		return(lapply(loci, collapse, allele=allele, ambig=ambiguity_code))
	} else {
		return(collapse(loci, allele=allele, ambig=ambiguity_code))
	}

}

