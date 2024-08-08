#' Clades collapsed to branches.
#' 
#' @description
#' Collapses specified clades to single branches.
#' 
#' @param tree an object of class `phylo` or `multiPhylo` or a character string with tree in newick format
#'   or a name of file containing the tree(s).
#' @param clades a data frame or matrix listing tip labels (1st column) and clade labels (2nd column).
#' @details The specified clades must be monophyletic in the tree(s), but may cover only a part of the tree(s).
#' @returns An object of class `phylo` of `multiPhylo` specifying a tree with a single tip per clade
#'   and terminal branch lengths equal to the mean root-to-tip distances in the clades.
#' @export

collapse_clades <- function(tree, clades) {
	find_common_anc <- function(tip, phy) {
		ifelse(length(tip) == 1, phy$edge[match(tip, phy$edge[,2]),1], ape::getMRCA(phy, tip))
	}
	find_clade_edges <- function(phy, from, to) {
		paths <- lapply(to, function(x) ape::nodepath(phy, from=from, to=x))
		paths <- lapply(paths, function(p) cbind(p[-length(p)], p[-1]))
		paths <- lapply(paths, function(p) apply(apply(p, 1, sort), 2, paste, collapse="-"))
		edges <- apply(apply(phy$edge, 1, sort), 2, paste, collapse="-")
		return(which(edges %in% unlist(paths)))
	}
	tree <- import_tree(tree)
	clades <- split(clades[,1], clades[,2])
	tips <- sapply(clades, "[", 1)
	drop <- setdiff(unlist(clades), tips)
	if (inherits(tree, "phylo")) {
		tree <- list(tree)
	}
	tree <- lapply(tree, function(x) ape::read.tree(text=ape::write.tree(ape::rotateConstr(phy=x, constraint=x$tip.label))))
	for (i in seq_along(tree)) {
		phy <- tree[[i]]	
		off <- lapply(clades, match, table=phy$tip.label)
		anc <- sapply(off, find_common_anc, phy=phy)
		edg <- lapply(seq_along(off), function(j) find_clade_edges(phy=phy, from=anc[j], to=off[[j]]))
		if (!is.null(phy$edge.length)) {
			brlen <- phy$edge.length[match(tips, phy$tip.label)]	
			brlen <- sapply(edg, function(jj) mean(phy$edge.length[jj])) - brlen
		}
		tips <- sapply(clades, "[", 1)
		tree[[i]] <- ape::drop.tip(phy, tip=drop)
		ord <- match(tips, tree[[i]]$tip.label)
		tree[[i]]$tip.label[ord] <- names(clades)
		terminal <- match(ord, tree[[i]]$edge[,2])
		if (!is.null(phy$edge.length)) {
			tree[[i]]$edge.length[terminal] <- tree[[i]]$edge.length[terminal] + brlen
		}
	}
	if (length(tree) == 1) {
		tree <- tree[[1]]
	}
	return(tree)
}

