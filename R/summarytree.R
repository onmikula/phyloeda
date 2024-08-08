#' Summary tree.
#' 
#' @description
#' Finds a summary tree for a posterior sample from a Bayesian phylogenetic inference.
#'
#' @param trees a posterior sample of rooted trees in a `multiPhylo` object or a list of `phylo` objects
#'   or a name of file with multiple trees in a format accepted by [read_tree].
#' @param type the type of summary tree, either `"MAP"` (maximum aposteriori, default),
#'   `"MCC"` (maximum clade credibility) or `"CON"` (majority consensus), see Details for explanantion.
#' @param file a name of the output file.
#' @param brlen logical, whether to estimate also branch lengths from the sample.
#' @param method string, name of the function to calculate branch lengths.
#' @param strict logical, whether to estimate branch lengths only from the trees, whose topology is strictly
#'   identical with that of a summary tree.
#' @param thin a proportion of mcmc steps to be retained or a subsampling frequency (if `thin > 1`).
#' @param digits the number of digits in PP values.
#' @details The maximum a posteriori tree is the most frequently sampled tree.
#'   The maximum clade credibility tree scores every sampled tree by the product of PPs of its constituent clades.
#'   The majority consensus tree adds the most supported clades unless they overlap (they are in conflict)
#'   with any of the clades that were already included. 
#' @returns An object of class `phylo` with tree whose node labels indicate posterior probabilities of the clades.
#' @export
                                                                                                                                                                                                                                                                                                                                  
summarytree <- function(trees, type=c("MAP","MCC","CON"), file, brlen=TRUE, method="mean", strict=FALSE, thin=1, digits=2) {
	if (is.character(trees) & length(trees) == 1) {
		trees <- read_tree(trees)
	}
	if (thin != 1) {
		trees <- thinning(trees, freq=thin)
	}
	tree_pp <- treePP(trees)
	summtrees <- setNames(vector("list", length(type)), type)
	for (i in seq_along(type)) {
		tree_fun <- match.fun(paste0("tree", toupper(type[i])))
		summtrees[[i]] <- tree_fun(tree_pp, digits=digits)
		if (isTRUE(brlen)) {
			summtrees[[i]] <- estimate_brlen(summtrees[[i]], trees, strict=strict, fun=method, thin=1)
		}		
	}
	if (length(summtrees) > 1) {
		class(summtrees) <- "multiPhylo"
	} else {
		summtrees <- summtrees[[1]]
	}
	if (!missing(file)) {
		ape::write.tree(summtrees, file)
	}
	return(summtrees)
}


#' @export                                                                                                                                                                                                                                                                                                                                   
treePP <- function(trees) {
	topocode <- function(phy, ref) {
		clades <- ape::prop.part(phy)
		labels <- attr(clades, "labels")
		clades <- cbind(diag(length(labels)), sapply(clades, function(x) as.numeric(ref %in% labels[x])))
		clades1 <- apply(-1 * (clades[,phy$edge[,2]] - 1), 2, paste, collapse="")
		clades2 <- apply(clades[,phy$edge[,2]], 2, paste, collapse="")
		biparts <- paste(sort(apply(cbind(clades1, clades2), 1, sort)[1,]), collapse="-")
		return(biparts)
	}
	nametips <- function(part, ref=attr(part, "labels")) lapply(lapply(part, function(x) ref[x]), sort)

	topolcodes <- unname(sapply(trees, topocode, ref=trees[[1]]$tip.label))
	topoltable <- sort(table(topolcodes), decreasing=TRUE)
	topologies <- trees[match(names(topoltable), topolcodes)]
	topologies <- unname(sapply(lapply(topologies, reorder_tree, tiporder=topologies[[1]]$tip.label), topology))
	treepp <- topoltable / length(trees)
	treedf <- data.frame(tree=topologies, pp=as.numeric(treepp))

	clades <- ape::prop.part(trees)
	cladepp <- attr(clades, "number") / length(trees)
	clades <- nametips(clades)
	if (length(treedf$tree) == 1) {
		credib <- list(nametips(ape::prop.part(ape::read.tree(text=treedf$tree))))
	} else {
		credib <- lapply(lapply(ape::read.tree(text=treedf$tree), ape::prop.part), nametips)
	}
	credib <- lapply(credib, match, table=clades)
	treedf$credib <- sapply(seq_along(credib), function(i) sum(log(cladepp[credib[[i]]])))
	clades <- sapply(clades, function(x) paste0("(", paste(x, collapse=","), ")"))
	cladedf <- data.frame(clade=clades, pp=cladepp)
	cladedf <- cladedf[order(cladedf$pp, decreasing=TRUE),]

	rooted <- ape::is.rooted(trees[[1]])

	result <- list(trees=treedf, clades=cladedf, rooted=rooted)
	return(result)
}


#' @export                                                                                                                                                                                                                                                                                                                                   
treeMAP <- function(pp, digits=2) {
	map <- ape::read.tree(text=pp$trees$tree[which.max(pp$trees$pp)])
	mapclades <- ape::prop.part(map)
	mapclades <- lapply(mapclades, function(x, labels) labels[x], labels=attr(mapclades, "labels"))
	mapclades <- lapply(mapclades, sort)
	mapclades <- sapply(mapclades, function(x) paste0("(", paste(x, collapse=","), ")"))
	map$node.label <- formatC(pp$clades$pp[match(mapclades, pp$clades$clade)], format="f", digits=digits)
	return(map)
}


#' @export                                                                                                                                                                                                                                                                                                                                   
treeMCC <- function(pp, digits=2) {
	mcc <- ape::read.tree(text=pp$trees$tree[which.max(pp$trees$cred)])
	mccclades <- ape::prop.part(mcc)
	mccclades <- lapply(mccclades, function(x, labels) labels[x], labels=attr(mccclades, "labels"))
	mccclades <- lapply(mccclades, sort)
	mccclades <- sapply(mccclades, function(x) paste0("(", paste(x, collapse=","), ")"))
	mcc$node.label <- formatC(pp$clades$pp[match(mccclades, pp$clades$clade)], format="f", digits=digits)
	return(mcc)
}


#' @export                                                                                                                                                                                                                                                                                                                                   
treeCON <- function(pp, digits=2) {
	clades <- strsplit(gsub("\\(|\\)", "", pp$clades$clade), ",")
	i <- 1
	con <- clades[i]
	while (i < length(clades)) {
		i <- i + 1
		included <- sapply(con, function(old, new) all(new %in% old), new=clades[[i]])
		includes <- sapply(con, function(old, new) all(old %in% new), new=clades[[i]])
		disjoint <- sapply(con, function(old, new) length(intersect(old, new)) == 0, new=clades[[i]])
		if (all(included | includes | disjoint)) {
			con <- append(con, clades[i])
		}
	}
	con <- con[order(sapply(con, length))]
	node.label <- formatC(pp$clades$pp[match(con, clades)], format="f", digits=digits)
	while (length(con) > 1) {
		nwk <- paste0("(", con[[1]][1], ",", con[[1]][2], ")", node.label[1])
		for (i in 2:length(con)) {
			ii <- match(con[[1]][1], con[[i]])
			if (!is.na(ii)) {
				con[[i]][ii] <- nwk
				con[[i]] <- setdiff(con[[i]], con[[1]])
			}
		}
		node.label <- node.label[-1]
		con <- con[-1]
	}
	nwk <- paste0("(", paste(con[[1]], collapse=","), ")", node.label[1], ";")
	con <- ape::read.tree(text=nwk)
	return(con)	
}



#' Branch lengths.
#' 
#' @description
#' Finds a summary tree for a posterior sample from a Bayesian phylogenetic inference.
#'
#' @param phy a `phylo` object with the summary tree (branch lengths are omitted, even if present).
#' @param trees a posterior sample of rooted trees in a `multiPhylo` object or a list of `phylo` objects
#'   or an object of class `bpp` or a name of file with the posterior sample.
#' @param strict logical, whether to estimate branch lengths only from the trees, whose topology is strictly
#'   identical with that of a summary tree.
#' @param fun character string, name of the function estimating branch lengths (default is `"mean"`).
#' @param thin a proportion of mcmc steps to be retained or a subsampling frequency (if `thin > 1`).
#' @returns An object of class `phylo` amended by branch lengths estimated from the provided sample of trees.
#' @export

estimate_brlen <- function(phy, trees, strict=FALSE, fun="mean", thin=1) {

	topocode <- function(phy, ref) {
		clades <- ape::prop.part(phy)
		labels <- attr(clades, "labels")
		clades <- cbind(diag(length(labels)), sapply(clades, function(x) as.numeric(ref %in% labels[x])))
		clades1 <- apply(-1 * (clades[,phy$edge[,2]] - 1), 2, paste, collapse="")
		clades2 <- apply(clades[,phy$edge[,2]], 2, paste, collapse="")
		biparts <- paste(sort(apply(cbind(clades1, clades2), 1, sort)[1,]), collapse="-")
		return(biparts)
	}
	bipartcodes <- function(phy, ref) {
		ntip <- ape::Ntip(phy)
		bipart <- character(nrow(phy$edge))
		for (i in seq_along(bipart)) {
			start <- phy$edge[i,1]
			fstep <- phy$edge[i,2]
			while(any(fstep > ntip)) {
				start <- c(start, fstep[fstep > ntip])
				fstep <- setdiff(as.numeric(phy$edge[phy$edge[,1] %in% fstep | phy$edge[,2] %in% fstep,]), start)
			}
			fstep <- as.numeric(ref %in% phy$tip.label[fstep])
			if (sum(fstep) > ntip / 2) {
				fstep <- -1 * (fstep - 1)
			}
			bipart[i] <- paste(fstep, collapse="")
		}
		return(bipart)
	}
	cladecodes <- function(phy, ref) {
		clades <- ape::prop.part(phy)
		labels <- attr(clades, "labels")
		clades <- cbind(diag(length(labels)), sapply(clades, function(x) as.numeric(ref %in% labels[x])))
		clades <- apply(clades, 2, paste, collapse="")
		return(clades)
	}
	
	fun <- match.fun(fun)
	if (thin != 1) {
		trees <- thinning(trees, freq=thin)
	}
	if (isTRUE(strict)) {
		trees <- trees[sapply(trees, topocode) == topocode(phy)]
		if (length(trees) == 0) {
			stop("No tree in the sample has the same topology as the target tree. Consider 'strict=FALSE' option.")
		}
	}

	rooted <- ape::is.rooted(phy)
	ultrametric <- ifelse(rooted, ape::is.ultrametric(trees[[1]], tol=0.001, option=1), FALSE)
	if (!rooted | (rooted & !ultrametric)) {
		target <- bipartcodes(phy, ref=phy$tip.label)
		sample <- lapply(trees, bipartcodes, ref=phy$tip.label)
		edgelengths <- sapply(seq_along(trees), function(i) trees[[i]]$edge.length[match(target, sample[[i]])])
		phy$edge.length <-  apply(edgelengths, 1, fun, na.rm=TRUE)
	} else {
		target <- cladecodes(phy, ref=phy$tip.label)
		sample <- lapply(trees, cladecodes, ref=phy$tip.label)
		nodeheights <- sapply(seq_along(trees), function(i) ape::node.depth.edgelength(trees[[i]])[match(target, sample[[i]])])
		nodeheights <-  apply(nodeheights, 1, fun, na.rm=TRUE)
		phy$edge.length <-  nodeheights[phy$edge[,2]] - nodeheights[phy$edge[,1]]
	}
	
	return(phy)
	
}

