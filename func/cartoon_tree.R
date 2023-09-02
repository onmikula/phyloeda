### FUNCTION
# cartoon_tree
### ARGUMENTS
# phy: tree object of class 'phylo'
# clades: data frame (or matrix) with tip labels in the first column and clade assignment in the second. The classification need not to be comprehensive.
# color: a named vector indicating colors of the cartoon triangles with the names matching clade labels. The color of ancestral and unclassified branches can be specified by an element named 'bg' (if not specified, black color is used).   
# alpha: opacity of the of the colors filling the cartoon triangles
# lwd: the width of the branches of the plotted phylogeny (~ edge.width of plot.phylo)
# cex: size of species label (if 'show.clade.label' is TRUE)
# show.clade.label: whether to show the clade labels along 
# type: the type of phylogeny (as in plot.phylo, but currently only 'phylogram' and 'cladogram' is implemented)
# use.edge.length: whether to use the edge lengths of the phylogeny to draw the branches (as in plot.phylo)
# direction: direction of the tree (as in plot.phylo, but currently only 'rightwards' is implemented)
# device: "quartz", "x11" (or "X11") or "pdf", if NULL, the objects are plotted into the current device
# file: name of pdf file if device == "pdf"
# width:	width of graphical device
# height:	height of graphical device
# mai: size of outer margins in inches, recycled if necessary
# extension: numeric vector of length 1 or 4 indicating percent by which the range of node coordinates is extended: c(bottom, left, top, right)
# return: whether to return coordinates of outer sides of cartoon triangles (x, mean y, min y, max y)
# expr: an expression whose evoluation adds a new element (e.g. axis) to the figure


cartoon_tree <- function(phy, clades, color, alpha=0.5, lwd=2, cex=1, show.clade.label=TRUE, type="phylogram", use.edge.length=TRUE, direction="rightwards", device, file, width=7, height=7, mai=0.02, extension=4, return=FALSE, expr) {
	opacity <- function(color, alpha) {
		hex <- c("0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F")
		num2hex <- function(x, hex) paste0(hex[(x-x%%16)/16+1], hex[x%%16+1])
		alpha <- as.integer(255 * rep(alpha, length(color)))
		rgb <- rbind(grDevices::col2rgb(color), alpha)
		return(paste0("#", apply(apply(rgb, 2, num2hex, hex=hex), 2, paste, collapse="")))
	}
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

	ntip <- ape::Ntip(phy)
	coord <- ape::plotPhyloCoor(phy, type=type, use.edge.length=use.edge.length, direction="rightwards", node.pos=NULL, tip.height=NULL)
	clades <- droplevels(clades[clades[,1] %in% phy$tip.label,])
	species <- split(clades[,1], clades[,2])
	tips <- lapply(species, function(tip, phy) match(tip, phy$tip.label), phy=phy)
	mrcas <- lapply(tips, find_ancestor, phy=phy)
	stem <- find_edges(anc=ntip + 1, off=setdiff(unlist(mrcas), seq(ntip)), phy)
	backgr <- find_edges(anc=ntip + 1, off=which(!phy$tip.label %in% clades[,1]), phy)
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

	if (isTRUE(show.clade.label) | isTRUE(return)) {
		spxy <- coord[unlist(ifelse(singletons, tips, mrcas)),c(1,2,2,2)]
		rownames(spxy) <- names(mrcas)
		colnames(spxy) <- c("x", "y", "miny", "maxy")
		for (i in which(!singletons)) {
			spxy[i,1] <- mean(coord[tips[[i]],1])
			spxy[i,3] <- min(coord[tips[[i]],2])
			spxy[i,4] <- max(coord[tips[[i]],2])
		}	
	}

	xlim <- range(coord[,1])
	ylim <- range(coord[,2])
	extension <- rep_len(extension / 100, 4)
	xlim <- xlim + c(-1, 1) * extension[c(2,4)] * diff(xlim)
	ylim <- ylim + c(-1, 1) * extension[c(1,3)] * diff(ylim)

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

	if (isTRUE(show.clade.label)) {
		text(spxy, rownames(spxy), cex=cex, pos=4)
	}
	
	
	if (isTRUE(device == "pdf")) {
		dev.off()
	}

	if (isTRUE(return)) {
		return(spxy)
	}

}
