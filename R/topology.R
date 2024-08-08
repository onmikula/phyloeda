#' Tree topology.
#' 
#' @description
#' Writes tree topology in Newick format.
#' 
#' @param phy an object of class `phylo` or `multiPhylo` or anything accepted by [read_tree].
#' @param tiporder a vector giving tip labels in the desired order or a `phylo` object from which
#'   the tip order is taken. Partial constraint on the tip order is accepted.
#' @returns A character string with tree topology (no branch lengths) in Newick format.
#' @export

topology <- function(tree, tiporder=NULL) {
	if (is.character(tree)) {
		tree <- read_tree(tree)
	}
	if (inherits(tree, "phylo")) {
		tree <- list(tree)
	}
	if (!is.null(tiporder)) {
		tree <- lapply(tree, reorder_tree, tiporder=tiporder)
	}
	tree <- lapply(tree, function(x) replace(replace(x, "edge.length", NULL), "node.label", NULL))
	tree <- unname(sapply(tree, ape::write.tree))
	return(tree)
}
