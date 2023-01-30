# partitioning of a rooted tree by the maximum cross-section branch length criterion
# phy - class 'hclust' or 'phylo'

max_crosssection <- function(phy) {
	sc <- function(k, h, b) {
		mean(b[nh[,2] > k & nh[,1] <= k])
	}
	phy <- as.phylo(phy)
	ntip <- ape::Ntip(phy)
	node.depth.edgelength(phy)
	nh <- matrix(ape::node.depth.edgelength(phy)[phy$edge], ape::Nedge(phy), 2)
	knots <- sort(setNames(nh[,2], phy$edge[,2])[!duplicated(nh[,2])])
	knots <- setNames(c(0, knots), c(ntip+1, names(knots)))
	k <- knots[which.max(sapply(knots, sc, h=nh, b=phy$edge.length))]
	anc <- phy$edge[nh[,2] > k & nh[,1] <= k, 2]
	single <- Filter(function(x) length(x) > 0, anc[anc <= ntip])
	otus <- lapply(setdiff(anc, single), function(a) ape::extract.clade(phy, a)$tip.label)
	otus <- c(as.list(phy$tip.label[single]), otus)
	otus <- data.frame(Tip=unlist(otus), Group=rep(seq_along(otus), sapply(otus, length)))
	return(otus)
}

