#' Write trees.
#' 
#' @description
#' Writes phylogenetic trees into a file written in newick or nexus format or into a newick text string.
#' 
#' @param phy an object of class `phylo` or `multiPhylo`.
#' @param file a character string with the file name.
#' @param format character string, either `"newick"` (default) or `"nexus"`
#' @param translate logical, relevant for the nexus format, whether to use numeric representation of tip labels
#' @param tree.names logical, relevant for the nexus format, whether to keep tree names
#' @returns If `file` is not supplied, it returns character vector with newick text strings.
#' @export

write_tree <- function(phy, file, format="newick", translate=TRUE, tree.names=TRUE) {

	if (inherits(phy, "phylo")) {
		phy <- list(phy)
	}
	format <- grep(format[1], c("newick", "nexus"), value=TRUE)[1]

	if (format == "newick") {
		 phy <- sapply(phy, ape::write.tree, file="")
	} else if (format == "nexus") {
		taxa <- phy[[1]]$tip.label
		if (isTRUE(tree.names)) {
			nam <- names(phy)
		} else {
			nam <- NULL
		}	
		nexus <- phy <- sapply(phy, ape::write.tree, file="")
		n <- length(taxa)
		header <- c("Begin taxa;", paste("\tDimensions ntax=", n, ";", sep=""), "\t\tTaxlabels", paste("\t\t", taxa, sep=""), "\t;", "End;")
		if (isTRUE(translate)) {
			transl <- paste(paste("\t\t", seq(n), sep=""), paste(taxa, ",", sep=""))
			transl[length(transl)] <- sub("\\,$", "", transl[length(transl)])
			transl <- c("\tTranslate", transl, ";")
			for (i in seq(n)) {
				pattern <- paste0(taxa[i], "[[:punct:]]{1}")
				m <- regexpr(pattern, nexus)
				attr(m, "match.length") <- attr(m, "match.length") - 1
				regmatches(nexus, m) <- i
			}
		} else {
			transl <- NULL
		}
		if (isTRUE(tree.names) & !is.null(nam)) {
			nexus <- paste("tree", nam, "=", nexus)
		}
		nexus <- c("#NEXUS", "\n", header, "\n", "Begin trees;", transl, nexus, "End;")
	}

	if (missing(file)) {
		return(phy)
	} else if (format == "newick") {
		writeLines(phy, con=file, sep="\n", useBytes=FALSE)
	} else if (format == "nexus") {
		writeLines(nexus, con=file, sep="\n", useBytes=FALSE)
	}
	
}

