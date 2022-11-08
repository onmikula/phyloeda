### FUNCTION
# paint_branches
### ARGUMENTS
# phy: tree object of class 'phylo'
# clades: Vector of clade labels (order must correspond to 'phy$tip.label') or data frame (or matrix) with tip labels in the first column and clade assignment in the second. The classification need not to be comprehensive.
# color: named vector of colors. Names correspond to either clade labels or tip labels (if coloration is unlinked from clades)
# stem: color of stem (unclassified) branches

paint_branches <- function(phy, clades, color, stem="black") {

	# functions
	getPaths <- function(from, to, phy) lapply(to, function(x) ape::nodepath(phy, from=from, to=x))
	getOffspring <- function(tip, phy) match(tip, phy$tip.label)
	getAncestor <- function(tip, phy) {
		ifelse(length(tip) == 1, phy$edge[phy$edge[,2] == match(tip, phy$tip.label), 1], ape::getMRCA(phy, tip))
	}
	makeEdges <- function(paths) {
		unique(do.call(rbind, lapply(paths, function(x) cbind(x[-length(x)], x[-1]))))
	}

	# defaults
	if (missing(clades)) {
		tips <- phy$tip.label
		clades <- rep_len(1, length(tips))
	} else if (!is.null(dim(clades))) {
		clades <- clades[match(phy$tip.label, clades[,1]),]
		tips <- clades[,1]
		clades <- clades[,2]
	} else if (is.vector(clades)) {
		tips <- phy$tip.label
		if (length(tips) != length(clades)) {
			stop("length of 'clades'must be equal to the number of tree tips")
		}
	}
	L <- sort(unique(clades))
	K <- length(L)
	N <- length(tips)
	if (missing(color)) {
		if (length(unique(clades)) == 1) {
			color <- setNames(stem, L)
		} else {
			palette <- c(Red="#e6194b", Green="#3cb44b", Yellow="#ffe119", Blue="#0082c8", Orange="#f58231", Purple="#911eb4", Cyan="#46f0f0", Magenta="#f032e6", Lime="#d2f53c", Pink="#fabebe", Teal="#008080", Lavender="#e6beff", Brown="#aa6e28", Beige="#fffac8", Maroon="#800000", Mint="#aaffc3", Olive="#808000", Coral="#ffd8b1", Navy="#000080", Grey="#808080")		
			color <- setNames(palette[seq(K)], L)	
		}
	} else if (is.null(names(color)) & length(color) >= K) {
		color <- setNames(color[seq(K)], L)
	} else if (any(names(color) %in% L)) {
		color <- setNames(color[L], L)
	} else if (any(names(color) %in% tips)) {
		color <- setNames(color[tips], tips)
	}
	if (is.numeric(color)) {
		color <- setNames(grDevices::palette()[color], names(color))
	}

	# coloring edges
	species <- split(tips, clades)
	off <- lapply(species, getOffspring, phy=phy)
	anc <- lapply(species, getAncestor, phy=phy)
	paths <- mapply(getPaths, anc, off, MoreArgs=list(phy=phy), SIMPLIFY=FALSE)	
	group <- lapply(paths, makeEdges)
	group <- lapply(group, function(x) match(x[,2], phy$edge[,2]))
	group <- data.frame(edge=unlist(group), otu=rep(names(group), sapply(group, length)), stringsAsFactors=FALSE)
	stems <- setdiff(seq(nrow(phy$edge)), group$edge)
	if (length(stems) > 0) {
		stems <- data.frame(edge=stems, otu="stem", stringsAsFactors=FALSE)
	} else {
		stems <- NULL
	}	
	group <- rbind(group, stems)
	group <- group[order(as.numeric(group$edge)),2]
	tipcolor <- group
	nodes <- sort(unique(c(phy$edge[group == "stem",1], phy$edge[group != "stem",2])))
	stems <- nodes[sapply(nodes, function(i) any(group[phy$edge[,1] == i] != group[phy$edge[,2] == i]))]
	while (length(stems) > 0) {
		group[(phy$edge[,1] %in% stems | phy$edge[,2] %in% stems)] <- "stem"
		nodes <- sort(unique(c(phy$edge[group == "stem",1], phy$edge[group != "stem",2])))
		stems <- nodes[sapply(nodes, function(i) any(group[phy$edge[,1] == i] != group[phy$edge[,2] == i]))]
	}
	group <- ifelse(phy$edge[,2] <= N, tipcolor, group)
	edge.color <- c(color, stem=stem)[group]
	return(edge.color)

}



