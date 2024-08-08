#' Branch colors.
#' 
#' @description
#' Assigns colors to branches depending which clade they belong to.
#' 
#' @param phy an object of class `phylo`.
#' @param info a data frame (or matrix) with tip labels (1st column) and clade assignment (2nd column).
#'   The classification need not to be comprehensive.
#' @param color named vector of colors, the names correspond to clade labels,
#'   or a matrix or data frame with colors (1st column) and clade labels (2nd column).
#' @param bg color of unclassified ("background") branches.
#' @param nested logical, whether to paint branches belonging to the nested clades (if `TRUE`)
#'   or retain the color of the larger clade (if `FALSE`).
#' @param returns A vector of colors assigned to tree branches, ordered as in `phy$edge`.
#' @export

paint_branches <- function(phy, info, color, bg="black", nested=TRUE) {
	nodepaths <- function(pairs, phy) {
		lapply(seq(nrow(pairs)), function(i) ape::nodepath(phy, from=pairs[i,1], to=pairs[i,2]))
	}
	makedges <- function(path) {
		apply(apply(cbind(path[-length(path)], path[-1]), 1, sort), 2, paste, collapse="-")
	}
	info <- info[match(intersect(phy$tip.label, info[,1]), info[,1]),,drop=FALSE]
	clades <- split(match(info[,1], phy$tip.label), info[,2])
	clades <- clades[order(sapply(clades, length), decreasing=TRUE)]
	if (!is.null(dim(color))) {
		color <- setNames(color[,1], color[,2])
	}
	color <- color[names(clades)]
	if (any(is.na(color))) {
		clades <- clades[!is.na(color)]
		color <- clades[!is.na(color)]
		warning("for some clades no color was specified")
	}
	edges <- apply(apply(phy$edge, 1, sort), 2, paste, collapse="-")
	edge.color <- rep(bg, length(edges))
	overlapping <- rep(FALSE, length(edges))
	paths <- lapply(lapply(lapply(clades, combn, m=2), t), nodepaths, phy=phy)
	for (i in seq_along(paths)) {
		paths[[i]] <- sort(unique(unlist(lapply(paths[[i]], makedges))))
	}
	for (i in seq_along(paths)) {
		ii <- match(paths[[i]], edges)
		colored <- !is.na(edge.color[ii])
		if (all(colored)) {
			if (isFALSE(nested)) {
				next
			}
		} else if (any(colored)) {
			overlapping[ii[colored]] <- TRUE
		}
		edge.color[ii] <- color[i]
	}
	edge.color[overlapping] <- bg
	return(edge.color)
}
