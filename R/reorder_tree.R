#' Reorder tree tips.
#' 
#' @description
#' Rotates clades to achieve specified ordering of tree tips (if possible, given the tree topology).
#'
#' @param phy an object of class `phylo` with the input tree.
#' @param tiporder a vector giving tip labels in the desired order or a `phylo` object from which
#'   the tip order is taken. Partial constraint on the tip order is accepted.
#' @returns Rotated tree as a `phylo` object.
#' @details The function uses iteratively [ape::rotate] and then exports the rotated `phy` object using
#'   `ape::read.tree(text=ape::write.tree(phy))`. If the input contains components `edge.label` or `edge.color`,
#'   they are retained and their elements reordered as necessary.
#' @export

reorder_tree <- function(phy, tiporder) {

	if (missing(tiporder)) {
		tiporder <- phy$tip.label[order(match(seq(ape::Ntip(phy)), phy$edge[,2]))]
	} else if (inherits(tiporder, "phylo")) {
		tiporder <- tiporder$tip.label[order(match(seq(ape::Ntip(tiporder)), tiporder$edge[,2]))]
	}
	tiporder <- intersect(tiporder, phy$tip.label)
	if (any(!phy$tip.label %in% tiporder)) {
		tiporder <- unique(c(tiporder, phy$tip.label[order(match(seq(ape::Ntip(phy)), phy$edge[,2]))]))
	}	
	edge.label <- phy$edge.label
	edge.color <- phy$edge.color
	edge.attr <- !is.null(edge.label) | !is.null(edge.color)
	if (isTRUE(edge.attr)) {
		edge_input <- phy$edge
	}
	
	clades <- c(as.list(phy$tip.label), lapply(ape::prop.part(phy), function(x) phy$tip.label[x]))
	constr <- lapply(clades, match, table=tiporder)
	nodes <- seq(ape::Nnode(phy)) + ape::Ntip(phy)
	for (i in nodes) {
		j <- phy$edge[phy$edge[,1] == i, 2]
		comb <- combn(length(j), 2, simplify=FALSE)
		for (k in seq_along(comb)) {
			jj <- comb[[k]]
			if (max(constr[j][[jj[1]]]) > min(constr[j][[jj[2]]])) {
				phy <- ape::rotate(phy, node=i, polytom=comb[[k]])
				j[comb[[k]]] <- rev(j[comb[[k]]])
			}
		}
	}

	if (isTRUE(edge.attr)) {
		ord <- base::match(apply(phy$edge, 1, paste, collapse="-"), apply(edge_input, 1, paste, collapse="-"))
	}
	phy$tip.label <- gsub(" ", "@", phy$tip.label)
	phy <- ape::read.tree(text=ape::write.tree(phy))
	phy$tip.label <- gsub("@", " ", phy$tip.label)
	if (isTRUE(edge.attr)) {
		if (!is.null(edge.label)) {
			phy$edge.label <- edge.label[ord]
		}
		if (!is.null(edge.color)) {
			phy$edge.color <- edge.color[ord]
		}
	}

	return(phy)

}



#' Reorder clades.
#' 
#' @description
#' Reorders clades in the tree, so they are plotted in the specified order.
#' 
#' @param phy an object of class `phylo` with the input tree.
#' @param clades a data frame listing tip labels (1st column) and clade labels (2nd column).
#'   The specified clades must be monophyletic in the tree, but may cover only a part of it.
#' @param cladeorder a vector giving clade labels in the desired order.
#' @returns Rotated tree as a `phylo` object.
#' @export

reorder_clades <- function(phy, clades, cladeorder) {
	clades <- clades[clades[,1] %in% phy$tip.label,]
	if (missing(cladeorder)) {
		cladeorder <- unique(clades[,2])
	} else {
		cladeorder <- intersect(unique(cladeorder), clades[,2])
	}
	tiporder <- unname(unlist(split(clades[,1], clades[,2])[cladeorder]))
	phy <- reorder_tree(phy, tiporder)
	return(phy)
}

