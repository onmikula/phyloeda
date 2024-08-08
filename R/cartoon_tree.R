#' Cartoon tree.
#' 
#' @description
#' Plots phylogenetic tree while collapsing defined clades into triangles.
#' 
#' @param phy an object of class `phylo`.
#' @param cartoon a data frame (or matrix) with tip labels (1st column) and clade assignment (2nd column).
#'   The classification need not to be comprehensive.
#' @param color named vector of colors, the names correspond to clade labels.
#'   Optionally, an element called `"bg"` can be included to define color of unclassified branches.
#' @param alpha numeric, opacity of triangle filling.
#' @param lwd width of the lines.
#' @param cex size of species labels.
#' @param show.sp.label logical, whether to show species labels.
#' @param type character, tree type, either `"phylogram"` (default) or `"cladogram"`.
#' @param use.edge.length logical, whether to use edge lengths. 
#' @param direction a character string specifyying direction of the tree.
#' @param device Either `"quartz"`, `"x11"` or `"pdf"` to indicate the type of graphical device.
#'   If `NULL`, the objects are plotted into the current device.
#' @param file a name of pdf file if `device == "pdf"`.
#' @param width width of graphical device in inches.
#' @param height height of graphical device in inches.
#' @param mai size of outer margins in inches, recycled if necessary.
#' @param xylim percents of node coordinate ranges to be added to x and y axes at their extremes.
#' @param return logical, whether to return coordinates of species labels.
#' @export

cartoon_tree <- function(phy, info, color, alpha=0.5, lwd=2, cex=1, show.sp.label=TRUE, type="phylogram", use.edge.length=TRUE, direction="rightwards", device, file, width=7, height=7, mai=0.02, xylim=4, return=FALSE) {

	find_ancestor <- function(tip, phy) ifelse(length(tip) == 1, phy$edge[match(tip, phy$edge[,2]),1], ape::getMRCA(phy, tip))
	find_edges <- function(anc, off, phy) {
		paths <- lapply(off, function(to) ape::nodepath(phy, from=anc, to=to))
		edges <- unique(do.call(rbind, lapply(paths, function(x) cbind(x[-length(x)], x[-1]))))
		edges <- paste(edges[,1], edges[,2], sep="-")
		return(match(edges, paste(phy$edge[,1], phy$edge[,2], sep="-")))
	}
	draw_phyloline <- function(xy, col, lwd)	 {
		lines(xy[c(1,1)], xy[c(3,4)], col=col, lwd=lwd)
		lines(xy[c(1,2)], xy[c(4,4)], col=col, lwd=lwd)
	}
	draw_cladoline <- function(xy, col, lwd)	 {
		lines(xy, col=col, lwd=lwd)
	}
	draw_line <- list(phylogram=draw_phyloline, cladogram=draw_cladoline)[[type]]
	draw_triangle <- function(anc, off, coord, col, bg, lwd) {
		xy <- rbind(coord[anc,], cbind(mean(coord[off,1]), range(coord[off,2])))
		polygon(xy, border=col, col=bg, lwd=lwd)
	}

	phy <- ape::read.tree(text=ape::write.tree(ape::rotateConstr(phy, constraint=phy$tip.label)))
	
	ntip <- ape::Ntip(phy)
	coord <- ape::plotPhyloCoor(phy, type=type, use.edge.length=use.edge.length, direction="rightwards", node.pos=NULL, tip.height=NULL)
	species <- split(info[,1], info[,2])
	tips <- lapply(species, function(tip, phy) match(tip, phy$tip.label), phy=phy)
	mrcas <- lapply(tips, find_ancestor, phy=phy)
	stem <- find_edges(anc=ntip + 1, off=setdiff(unlist(mrcas), seq(ntip)), phy)
	backgr <- find_edges(anc=ntip + 1, off=which(!phy$tip.label %in% info[,1]), phy)
	singletons <- sapply(species, length) == 1
	for (i in which(!singletons)) {
		coord[mrcas[[i]],2] <- mean(range(coord[tips[[i]],2]))
	}
	
	if (missing(color)) {
		color <- setNames(rep("#000000", length(species) + 1), c(names(species), "bg"))
	}
	if (!"bg" %in% names(color)) {
		color <- c(color, bg="#000000")
	}
	missingsp <- setdiff(names(species), names(color))
	color <- c(color, setNames(rep(color["bg"], length(missingsp)), missingsp))
	bgcolor <- setNames(opacity(color, alpha), names(color))

	if (isTRUE(show.sp.label) | isTRUE(return)) {
		spxy <- coord[unlist(ifelse(singletons, tips, mrcas)),]
		rownames(spxy) <- names(mrcas)
		for (i in which(!singletons)) {
			spxy[i,1] <- mean(coord[tips[[i]],1])
		}	
	}

	xlim <- range(coord[,1])
	ylim <- range(coord[,2])
	xylim <- rep_len(xylim / 100, 4)
	xlim <- xlim + c(-1, 1) * xylim[c(2,4)] * diff(xlim)
	ylim <- ylim + c(-1, 1) * xylim[c(1,3)] * diff(ylim)

	if (missing(device)) {
		device <- ifelse(.Platform$OS.type == "unix", "quartz", "x11")
	}
	if (isTRUE(device == "pdf")) {
		pdf(ifelse(missing(file), "tree.pdf", file), width=width, height=height)
	} else if (!is.null(device)) {
		match.fun(tolower(device))(width=width, height=height)	
	}

	par(mai=rep_len(mai, 4), xaxs="i", yaxs="i")
	plot(coord, type="n", axes=FALSE, ann=FALSE, xlim=xlim, ylim=ylim)
	for (i in stem) {
		draw_line(coord[phy$edge[i,],], col=color["bg"], lwd=lwd)
	}
	for (i in backgr) {
		draw_line(coord[phy$edge[i,],], col=color["bg"], lwd=lwd)
	}
	for (i in which(singletons)) {
		branch <- bottom <- coord[phy$edge[match(tips[[i]], phy$edge[,2]),],]
		bottom[2] <- bottom[2] - 0.2 * diff(bottom[,1])
		top <- rbind(bottom[2,], branch[2,])
		draw_line(bottom, col=color["bg"], lwd=lwd)		
		draw_line(top, col=color[names(species)[i]], lwd=lwd)		
	}
	for (i in which(!singletons)) {
		sp <- names(species)[i]
		draw_triangle(anc=mrcas[[i]], off=tips[[i]], coord=coord, col=color[sp], bg=bgcolor[sp], lwd=lwd)
	}
	if (isTRUE(show.sp.label)) {
		text(spxy, rownames(spxy), cex=cex, pos=4)
	}
	
	if (isTRUE(device == "pdf")) {
		invisible(dev.off())
	}

	if (isTRUE(return)) {
		return(spxy)
	}

}
