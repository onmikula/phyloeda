remove_repeats <- function(loci, minn=5, kmer=5, delete=FALSE) {

	rm_rep <- function(dna, minn, kmer, delete) {

		dimers <- c("CA","TA","GA","AC","TC","GC","AT","CT","GT","AG","CG","TG")
		if (is.character(dna) & length(dna) == 1) {
			dna <- toupper(ape::read.dna(dna, format="fasta", as.matrix=TRUE, as.character=TRUE))
		} else if (inherits(dna, "DNAStringSet")) {
			dna <- toupper(as.character(as.matrix(ape::as.DNAbin(dna))))
		} else if (inherits(dna, "DNAbin")) {
			dna <- toupper(as.character(as.matrix(dna)))
		} else {
			dna <- toupper(dna)
		}
	
		dna[!dna %in% c("A","C","T","G")] <- "-"
		char <- apply(dna, 1, paste, collapse="")
		n <- length(char)
	
		dim <- setNames(lapply(dimers, gregexpr, text=char), dimers)
		dim <- lapply(dim, function(d) lapply(d, as.numeric))
		for (i in seq_along(dim)) {
			for (j in seq(n)) {
				dij <- dim[[i]][[j]]
				if (dij[1] != -1) {
					dij <- cbind(dij[-length(dij)], dij[-1]-1)
					dij <- dij[(dij[,2] - dij[,1] + 1) <= kmer,,drop=FALSE]
					nr <- nrow(dij)
					if (nr > 0) {
						dim[[i]][[j]] <- unname(unique(sapply(seq(nr), function(k) substr(char[j], dij[k,1], dij[k,2]))))
					} else {
						dim[[i]][[j]] <- list()
					}		
				} else {
					dim[[i]][[j]] <- list()
				}
			}
		}
		motifs <- sort(unique(unlist(dim)))
		motifs <- motifs[!grepl("-", motifs)]

		mot <- setNames(lapply(paste0("(", motifs, "){", minn, ",}"), gregexpr, text=char), motifs)
		for (i in seq_along(mot)) {
			for (j in seq(n)) {
				if (mot[[i]][[j]][1] != -1) {
					r <- list(mot[[i]][[j]], mot[[i]][[j]] + attr(mot[[i]][[j]],"match.length") - 1)
					r <- Map(seq, r[[1]], r[[2]])
					dna[j,sort(unique(unlist(r)))] <- "-"
				}
			}
		}
	
		if (isTRUE(delete)) {
			dna <- dna[,!apply(dna == "-", 2, all), drop=FALSE]
		}

		return(dna)
	
	}

	if (is.list(loci)) {
		bedtable <- attr(loci, "bedtable")
		loci <- lapply(loci, rm_rep, minn=minn, kmer=kmer, delete=delete)
		attr(loci, "bedtable") <- bedtable
	} else {
		rm_rep(dna=loci, minn=minn, kmer=kmer, delete=delete) 
	}

	return(loci)
	
}
