# X - input file with phylogenetic trees or DNA alignment
#	trees: object of class "phylo", "multiPhylo" or a list of "phylo" objects
#	DNA: character matrix or object of class "DNAbin" 
# file - name of output file

writeNexus <- function(X, file, translate=TRUE, tree.names=TRUE, part=NULL, interleave=FALSE, datatype="DNA", gap="-", missing="?") {
	cls <- class(X)
	if (cls[1] == "list") {
		cls <- class(X[[1]])
	}
	if (cls[1] %in% c("phylo", "multiPhylo")) {
		if (class(X) == "phylo") {
			taxa <- X$tip.label
			trees <- write.tree(X, file="")
		} else {
			taxa <- X[[1]]$tip.label
			trees <- unlist(lapply(X, write.tree, file=""))
		}
		n <- length(taxa)
		header <- c("Begin taxa;", paste("\tDimensions ntax=", n, ";", sep=""), "\t\tTaxlabels", paste("\t\t", taxa, sep=""), "\t;", "End;")
		if (isTRUE(translate)) {
			transl <- paste(paste("\t\t", seq(n), sep=""), paste(taxa, ",", sep=""))
			transl[length(transl)] <- sub("\\,$", "", transl[length(transl)])
			transl <- c("\tTranslate", transl, ";")
			for (i in seq(n)) {
				pattern <- paste0(taxa[i], "[[:punct:]]{1}")
				m <- regexpr(pattern, trees)
				attr(m, "match.length") <- attr(m, "match.length") - 1
				regmatches(trees, m) <- i
			}
		} else {
			transl <- NULL
		}
		if (isTRUE(tree.names)) {
			if(class(X) %in% c("list","multiPhylo") & !is.null(names(X))) {
				nam <- names(X)
			} else {
				nam <- paste("no", seq_along(trees), sep="_")
			}
			trees <- paste("tree", nam, "=", trees)
		}
		nexus <- c("#NEXUS", "\n", header, "\n", "Begin trees;", transl, trees, "End;")
	} else if (cls[1] %in% c("matrix", "DNAbin")) {
		seqcollapse <- function(M, collapse="") unname(apply(M, 1, paste, collapse=collapse))
		if (is.null(part)) {
			part <- attr(X, "part")
		}
		if (cls[1] == "DNAbin") {
			X <- as.matrix(X)
		}
		dims <- dim(X)
		taxa <- rownames(X)
		format <- paste("FORMAT", paste0("DATATYPE=", datatype), paste0("GAP=", gap), paste0("MISSING=", missing, ";"))
		header <- c("#NEXUS\n", "BEGIN DATA;",
			paste("DIMENSIONS", "  ", "NTAX=", dims[1], " ", "NCHAR=", dims[2], ";", sep=""),
			format, "MATRIX\n")
		if (isTRUE(interleave) & dims[2] > 80) {
			header[4] <- "FORMAT DATATYPE=DNA INTERLEAVE GAP=~ MISSING=-;"
			intervals <- findInterval(seq(dims[2]), seq(1, dims[2], 80))
			sequences <- lapply(seq(max(intervals)), function(i) cbind(taxa, seqcollapse(X[,intervals == i,drop=FALSE])))
			sequences <- lapply(sequences, seqcollapse, collapse=" ")
			sequences <- unlist(Map(c, sequences, ""))
		} else {
			sequences <- c(taxa, seqcollapse(X))
			sequences <- sequences[rep(seq(dims[1]), each=2) + (seq(2 * dims[1]) %% 2 == 0) * dims[1]]			
		}
		nexus <- c(header, sequences, ";\n", "END;")
		if (!is.null(part)) {
			p <- apply(part, 1, function(x) paste(x, collapse=" - "))
			p <- paste(names(p), p, sep=" = ")
			p <- paste(paste(paste("\t", "CHARSET", sep=""), p), ";", sep="")
			p <- c("BEGIN assumptions;", p, "END;")
			nexus <- c(nexus, "\n", p)
		}
	}
	if (missing(file)) {
		return(nexus)
	} else {
		writeLines(nexus, con=file, sep="\n", useBytes=FALSE)
	}
}

