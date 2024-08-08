#' Read trees.
#' 
#' @description
#' Reads phylogenetic trees from a newick text string or a file written in newick or nexus format.
#' 
#' @param tree a character string with a tree written in newick format or a file name,
#'   or an object of class `phylo` or `multiPhylo`.
#' @param thin a proportion of mcmc steps to be retained or a subsampling frequency (if `thin > 1`).
#' @returns An object of class `phylo` or `multiPhylo`.
#' @export

read_tree <- function(tree, thin=1) {
	if (inherits(tree, "phylo") | inherits(tree, "multiPhylo")) {
		phy <- tree
	} else if (all(sapply(tree, inherits, what="phylo"))) {
		phy <- tree
	} else {
		newick <- all(sapply(tree, function(x) grepl("^\\(.+;$", x)))
		if (isTRUE(newick)) {
			if (length(tree) > 1) {
				phy <- sapply(tree, function(text) ape::read.tree(text=text))
			} else {
				phy <- ape::read.tree(text=tree)
			}
		} else {
			first <- gsub("^[[:blank:]]+|[[:blank:]]+$", "", base::readLines(tree, n=100))
			first <- first[nchar(first) > 0]
			newick <- grepl("^\\(.+;$", first[1])
			nexus <- any(grepl("nexus", first, ignore.case=TRUE))
			if (isTRUE(newick)) {
				phy <- ape::read.tree(file=tree)
			} else if (isTRUE(nexus)) {
				phy <- ape::read.nexus(file=tree)
			}
		}
	}
	if (all(sapply(phy, inherits, what="phylo"))) {
		class(phy) <- "multiPhylo"
	}
	if (inherits(phy, "multiPhylo") & thin != 1) {
		phy <- thinning(phy, freq=thin)
	}
	return(phy)
}

