#' Clades in the tree.
#' 
#' @description
#' Plots phylogenetic tree while highlighting defined clades.
#' 
#' @param phy an object of class `phylo`.
#' @param clades vector of clade labels (order must correspond to `phy$tip.label`) or a data frame (or matrix)
#'   with tip labels (1st column) and clade assignment (2nd column). The classification need not to be comprehensive.
#' @param color named vector of colors. The names correspond to either clade labels or tip labels, the latter is
#'   relevant especially if coloration is unlinked from clade assignment. It can be also a matrix or data frame 
#'   with labels (1st column) and colors (2nd column).
#' @param bg color of background (unclassified) branches.
#' @param type character, tree type, either `"phylogram"` (default) or `"cladogram"`.
#' @param edge.color branch colors, vector of the same length as `nrow(phy$edge)`, if specified,
#'   it supersedes `color` of branches, while clade lines and names are not affected.
#' @param edge.width numeric, branch width.
#' @param show.clade.lines logical, whether to show clade lines.
#' @param show.clade.label logical, whether to show clade names.
#' @param show.tip.label logical, whether to display tip labels.
#' @param show.node.label logical, whether to display node labels.
#' @param show.edge.label logical, whether to display edge labels.
#' @param tip_labels vector of tip labels, in order corresponding to `phy$tip.label`.
#' @param clade_labels a named vector of clade labels, the names are clade names;
#'   useful to modify clade names for display (e.g. include spaces, delete underscores).
#' @param node_labels vector of node labels, in order corresponding to node numbers in `phy$edge`.
#' @param edge_labels vector of edge labels, in order corresponding to `phy$edge`.
#' @param root_label a character string added as a root label to `node_labels`.
#' @param ancestral logical, whether to show `node_labels` and `edge_labels` only for clade ancestors.
#' @param clade.order logical (whether to order the clades alphabetically) or character
#'   (clade labels in the required order).
#' @param show.scale logical, whether to display scale bar as x-axis.
#' @param device either `"quartz"`, `"x11"` or `"pdf"` to indicate the type of graphical device.
#'   If `NULL`, the objects are plotted into the current device.
#' @param width width of graphical device in inches.
#' @param height height of graphical device in inches.
#' @param file a name of pdf file if `device == "pdf"`.
#' @param outer sizes of outer margins in inches, recycled if necessary.
#' @param inner sizes of inner margins as proportions of plot dimensions, recycled if necessary.
#' @param clade.lin.pos position of clade lines, negative values ~ shift to the left.
#' @param clade.lin.over overhang of clade lines as a proportion of single character height.
#' @param clade.lab.col color of clade names.
#' @param clade.lab.cex size of clade names.
#' @param clade.lab.font font of clade names.
#' @param clade.lab.las orientation of clades names relative to the axis (default perpendicular).
#' @param tip.lab.pos gap between tree tips and tip labels as a proportion of the plot width (without outer margins).
#' @param tip.lab.cex size of tip labels.
#' @param tip.lab.font font of tip labels.
#' @param node.lab.pos node label position as a proportion of the supporting branches.
#' @param node.lab.cex size of node labels.
#' @param node.lab.font font of node labels.
#' @param node.lab.col color of node labels.
#' @param node.lab.border color of node label point border.
#' @param node.lab.bg color of node label point background.
#' @param node.lab.lwd width of node label point border.
#' @param node.lab.asp shape of node label, aspect ratio of an ellipse, default=1 means cricle.
#' @param node.lab.adj adjustment of text position in node label (`adj` argument of [graphics::text]).
#' @param edge.lab.pos edge label position as a proportion of the branch from its start (= root side).
#' @param edge.lab.cex size of edge label point.
#' @param edge.lab.font font of edge labels.
#' @param edge.lab.col color of edge labels.
#' @param edge.lab.border color of edge label point border.
#' @param edge.lab.bg color of edge label point background.
#' @param edge.lab.lwd width of edge label point border.
#' @param edge.lab.asp shape of edge label, aspect ratio of an ellipse, default=1 means cricle.
#' @param edge.lab.adj adjustment of text position in edge label (`adj` argument of [graphics::text]).
#' @param scale.cex size of scale labels.
#' @export


plot_clades <- function(phy, clades, color, bg="black", type="phylogram", edge.color=NULL, edge.width=3, show.clade.lines=TRUE, show.clade.label=TRUE, show.tip.label=TRUE, show.node.label=FALSE, show.edge.label=FALSE, clade_labels=NULL, tip_labels=NULL, node_labels=NULL, edge_labels=NULL, root_label=NULL, ancestral=TRUE, clade.order=FALSE, show.scale=TRUE, device, width=7, height=7, file=NULL, outer=NULL, inner=NULL, clade.lin.pos=0, clade.lin.over=0.4, clade.lab.col=NULL, clade.lab.cex=1.25, clade.lab.font=1, clade.lab.las=2, tip.lab.pos=0.005, tip.lab.cex=1, tip.lab.font=3, node.lab.pos=0, node.lab.cex=1, node.lab.font=1, node.lab.col="black", node.lab.border="black", node.lab.bg="white", node.lab.lwd=1, node.lab.asp=1, node.lab.adj=c(0.5,0.5), edge.lab.pos=0.5, edge.lab.cex=1, edge.lab.font=1, edge.lab.col="black", edge.lab.border="black", edge.lab.bg="white", edge.lab.lwd=1, edge.lab.asp=1, edge.lab.adj=c(0.5,0.5), scale.cex=1, scale.pos=c(0.5,0)) {

	drawTree <- list(phylogram=drawPhylogram, cladogram=drawCladogram)[[type]]

	phy <- ape::read.tree(text=ape::write.tree(ape::rotateConstr(phy, constraint=phy$tip.label)))

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
			stop("length of 'clades' must be equal to the number of tree tips")
		}
	}
	L <- sort(unique(clades))
	K <- length(L)
	N <- length(tips)
	ntip <- length(setdiff(phy$edge[,2], phy$edge[,1]))
	if (missing(color)) {
		palette <- c(Red="#e6194b", Green="#3cb44b", Yellow="#ffe119", Blue="#0082c8", Orange="#f58231", Purple="#911eb4", Cyan="#46f0f0", Magenta="#f032e6", Lime="#d2f53c", Pink="#fabebe", Teal="#008080", Lavender="#e6beff", Brown="#aa6e28", Beige="#fffac8", Maroon="#800000", Mint="#aaffc3", Olive="#808000", Coral="#ffd8b1", Navy="#000080", Grey="#808080")
		color <- setNames(palette[seq(K)], L)
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
	if (!is.null(edge.color)) {
		edge.col <- edge.color
	} else {
		edge.col <- phyloeda::paint_branches(phy, info=cbind(tips, clades), color=color, bg=bg, nested=TRUE)
	}
	backgr <- edge.col == bg

	# plotting
	if (missing(device)) {
		device <- ifelse(.Platform$OS.type == "unix", "quartz", "x11")
	}
	if (isTRUE(device == "pdf") & missing(file)) {
		file <- "tree_clades.pdf"
	}
	if (isTRUE(device == "pdf")) {
		pdf(file, width=width, height=height)
	} else if (!is.null(device)) {
		match.fun(tolower(device))(width=width, height=height)
	}

# margins
	missmai <- is.null(outer)
	if (isTRUE(show.clade.label)) {
		strsize <- ifelse(clade.lab.las %in% c(1,2), strwidth, strheight)
		strw <- max(strsize(s=clades, units="inches", cex=1, font=clade.lab.font)) * clade.lab.cex
		onecharw <- mean(strw / nchar(clades))
	}
	if (isTRUE(missmai)) {
		mai <- c(0.02, 0.02, 0.02, 0.02)
		if (isTRUE(show.scale)) {
			onecharh <- strheight(s="O", units="inches", cex=1) * scale.cex
			mai[1] <- max(c(0.02, 3.5 * scale.pos[1] * 0.2 + 0.25 * onecharh))
		}
		if (isTRUE(show.clade.lines)) {
			mai[4] <- max(c(0.02, mai[4] + clade.lin.pos * 0.2))
		}
		if (isTRUE(show.clade.label)) {
			mai[4] <- mai[4] + (1.0 * onecharw + strw + 0.5 * onecharw)
		}
	} else {
		mai <- rep_len(outer, 4)
	}

# tip labels
	if (!is.null(tip_labels)) {
		tip_labels <- tip_labels[order(match(seq_along(phy$tip.label), phy$edge[,2]))]
	} else if (isTRUE(show.tip.label)) {
		tip_labels <- phy$tip.label[order(match(seq_along(phy$tip.label), phy$edge[,2]))]
	} else {
		tip_labels <- rep(NA, length(phy$tip.label))
	}
	labtips <- !is.na(tip_labels)
	ntiplab <- sum(labtips)
	
# node labels
	if (!is.null(node_labels)) {
		node_labels <- node_labels[order(c(0, match((2:phy$Nnode)+ntip, phy$edge[,2])))]
	} else if (isTRUE(show.node.label) & !is.null(phy$node.label)) {
		node_labels <- phy$node.label[order(c(0, match((2:phy$Nnode)+ntip, phy$edge[,2])))]
	} else {
		node_labels <- rep(NA, phy$Nnode)
	}
	if (!is.null(root_label)) {
		node_labels[1] <- root_label
	}
	labnodes <- !is.na(node_labels)
	if (isTRUE(ancestral)) {
		nodeord <- match((2:phy$Nnode)+ntip, phy$edge[,2])
		labnodes[2:phy$Nnode] <- ifelse(backgr[nodeord], labnodes[2:phy$Nnode], FALSE)
	}
	nnodelab <- sum(labnodes)
	node_labels <- node_labels[labnodes]
	node.lab.asp <- rep_len(node.lab.asp, nnodelab)
	node.lab.cex <- rep_len(node.lab.cex, nnodelab)
	node.lab.font <- rep_len(node.lab.font, nnodelab)
	node.lab.pos <- rep_len(node.lab.pos, nnodelab)

# edge labels
	if (isTRUE(show.edge.label) & is.null(edge_labels) & !is.null(phy$node.label)) {
		edge_labels <- phy$node.label[match(phy$edge[,2], seq(phy$Nnode) + ntip)]
	} else {
		edge_labels <- rep(NA, nrow(phy$edge))
	}
	labedges <- !is.na(edge_labels)
	if (isTRUE(ancestral)) {
		labedges <- ifelse(backgr, labedges, FALSE)
	}
	nedgelab <- sum(labedges)
	edge_labels <- edge_labels[labedges]
	edge.lab.asp <- rep_len(edge.lab.asp, nedgelab)
	edge.lab.cex <- rep_len(edge.lab.cex, nedgelab)
	edge.lab.font <- rep_len(edge.lab.font, nedgelab)
	edge.lab.pos <- rep_len(edge.lab.pos, nedgelab)

# layout - draft
	orig.edge.length <- phy$edge.length
	maxdist <- max(ape::dist.nodes(phy)[ntip+1,])
	phy$edge.length <- phy$edge.length / maxdist
	xy <- ape::plotPhyloCoor(phy, type=type, direction="rightwards", use.edge.length=TRUE)
	xy[,2] <- (xy[,2] - 1) / (ntip - 1)

# plot
	def.par <- par(no.readonly=TRUE)
	par(mai=mai, cex=1)
	plot(0.5, 0.5, xlim=c(0, 1), ylim=c(0, 1), xaxs="i", yaxs="i", ann=FALSE, axes=FALSE, type="n")

#abline(h=c(0,1), v=c(0,1))

# inner margins
	missinner <- is.null(inner)
	if (isTRUE(missinner)) {
		din <- par("din")
		ww <- din[1] - sum(mai[c(2,4)])
		hh <- din[2] - sum(mai[c(1,3)])
		prop <- 0.075 / min(c(ww, hh))
		inner <- rep(prop, 4)		
	} else {
		inner <- rep_len(inner, 4)
	}

# node & edge labels
	maxprop <- 0.05
	if (any(labnodes)) {
		nodestrw <- max(strwidth(node_labels, units="user", cex=1, font=node.lab.font) * node.lab.cex)
		if (nodestrw > maxprop) {
			node.lab.cex <- maxprop * node.lab.cex / nodestrw
		}
		node.lab.cex <- edge.lab.cex <- min(c(node.lab.cex, edge.lab.cex))
		nodestrw <- strwidth(node_labels, units="user", cex=1, font=node.lab.font) * node.lab.cex
		nchnode <- nchar(node_labels)
		node.lab.size <- rep(max((nchnode + 0.5) * mean(nodestrw / nchnode)), length(nodestrw))
	} else {
		node.lab.size <- 0
	}
	if (any(labedges)) {
		edgestrw <- max(strwidth(edge_labels, units="user", cex=1, font=edge.lab.font) * edge.lab.cex)
		if (edgestrw > maxprop) {
			edge.lab.cex <- maxprop * edge.lab.cex / edgestrw
		}
		edgestrw <- strwidth(edge_labels, units="user", cex=1, font=edge.lab.font) * edge.lab.cex
		nchedge <- nchar(edge_labels)
		edge.lab.size <- rep(max((nchedge + 0.5) * mean(edgestrw / nchedge)), length(edgestrw))
	} else {
		edge.lab.size <- 0
	}
	if (isTRUE(missinner) & (any(labnodes) | any(labedges))) {	
		if (any(labnodes)) {
			labnodexy <- xy[-seq(ntip),][labnodes,,drop=FALSE]
			if (any(node.lab.pos != 0)) {
				ord <- match(which(labnodes) + ntip, phy$edge[,2])
				edge_lengths <- xy[phy$edge[ord,2],1] - xy[phy$edge[ord,1],1]
				nonroot <- !is.na(edge_lengths)
				labnodexy[nonroot,1] <- labnodexy[nonroot,1] - (1 - node.lab.pos[nonroot]) * edge_lengths[nonroot]
			}
			below <- max(-labnodexy[,2] + 0.5 * node.lab.size)
			left <- max(-labnodexy[,1] + 0.5 * node.lab.size)
			above <- max(labnodexy[,2] + 0.5 * node.lab.size - 1)
		}
		if (any(labedges)) {
			labedgexy <- xy[phy$edge[,2],][labedges,,drop=FALSE]
			if (any(edge.lab.pos != 0)) {
				edge_lengths <- phy$edge.length[labedges]
				labedgexy[,1] <- labedgexy[,1] - (1 - edge.lab.pos) * edge_lengths
			}
			below <- max(c(below, -labedgexy[,2] + 0.5 * edge.lab.size))
			left <- max(c(left, -labedgexy[,1] + 0.5 * edge.lab.size))
			above <- max(c(above, labedgexy[,2] + 0.5 * edge.lab.size - 1))
		}
		inner[1] <- inner[1] + below * (below > 0)
		inner[2] <- inner[2] + left * (left > 0)
		inner[3] <- inner[3] + above * (above > 0)
	}
	tipstrh <- max(strheight(s=tip_labels[labtips], units="user", cex=1, font=tip.lab.font)) * tip.lab.cex
	yfactor <- (1 - inner[3] - inner[1])
	ycex <- 1.25 * tipstrh / (xy[2,2] * yfactor)
	if (ycex > 1) {
		tip.lab.cex <- tip.lab.cex / ycex
	}
	tipstrw <- strwidth(s=tip_labels[labtips], units="user", cex=1, font=tip.lab.font) * tip.lab.cex
	inner[4] <- inner[4] + max((xy[seq(ntip),1][labtips] + tipstrw + tip.lab.pos) - 1)
	xfactor <- (1 - inner[4] - inner[2])
	yfactor <- (1 - inner[3] - inner[1])
	xy[,1] <-  xy[,1] * xfactor + inner[2]
	xy[,2] <-  xy[,2] * yfactor + inner[1]

# tree branches
	drawTree(phy, xy, lwd=edge.width, col=edge.col)

# tip labels
	if (any(labtips)) {
		tipxy <- xy[seq(ntip),][labtips,,drop=FALSE]
		tipxy[,1] <- tipxy[,1] + tip.lab.pos
		text(tipxy, label=tip_labels[labtips], cex=tip.lab.cex, font=tip.lab.font, adj=c(0, 0.5))
	}

# node labels
	if (any(labnodes)) {
		nodexy <- xy[-seq(ntip),][labnodes,,drop=FALSE]
		node.lab.border <- rep_len(node.lab.border, nnodelab)
		node.lab.bg <- rep_len(node.lab.bg, nnodelab)
		node.lab.lwd <- rep_len(node.lab.lwd, nnodelab)
		node.lab.col <- rep_len(node.lab.col, nnodelab)
		if (!is.list(node.lab.adj)) {
			node.lab.adj <- list(node.lab.adj)
		}
		node.lab.adj <- rep_len(node.lab.adj, nnodelab)
		if (any(node.lab.pos != 0)) {
			ord <- match(which(labnodes) + ntip, phy$edge[,2])
			edge_lengths <- xy[phy$edge[ord,2],1] - xy[phy$edge[ord,1],1]
			nonroot <- !is.na(edge_lengths)
			nodexy[nonroot,1] <- nodexy[nonroot,1] - (1 - node.lab.pos[nonroot]) * edge_lengths[nonroot]
		}
		node.lab.asp <- node.lab.asp * (width - sum(mai[c(2,4)])) / (height - sum(mai[c(1,3)]))
		for (i in seq(nnodelab)) {
			drawEllipse(x=nodexy[i,1], y=nodexy[i,2], size=node.lab.size[i], w=1, h=node.lab.asp[i], r=0, res=200, col=node.lab.bg[i], border=node.lab.border[i], lwd=node.lab.lwd[i])		
			text(nodexy[i,,drop=FALSE], label=node_labels[i], cex=node.lab.cex[i], font=node.lab.font[i], adj=node.lab.adj[[i]], col=node.lab.col[i])
		}
	}

# edge labels
	if (any(labedges)) {
		edgexy <- xy[phy$edge[labedges,2],]
		edge.lab.col <- rep_len(edge.lab.col, nedgelab)
		edge.lab.border <- rep_len(edge.lab.border, nedgelab)
		edge.lab.bg <- rep_len(edge.lab.bg, nedgelab)
		edge.lab.lwd <- rep_len(edge.lab.lwd, nedgelab)
		if (!is.list(edge.lab.adj)) edge.lab.adj <- list(edge.lab.adj)
		edge.lab.adj <- rep_len(edge.lab.adj, nedgelab)
		if (any(edge.lab.pos != 0)) {
			edge_lengths <- phy$edge.length[labedges]
			edgexy[,1] <- edgexy[,1] - (1 - edge.lab.pos) * edge_lengths
		}
		edge.lab.asp <- edge.lab.asp * (width - sum(mai[c(2,4)])) / (height - sum(mai[c(1,3)]))
		for (i in seq(nedgelab)) {
			drawEllipse(x=edgexy[i,1], y=edgexy[i,2], size=edge.lab.size[i], w=1, h=edge.lab.asp[i], r=0, res=200, col=edge.lab.bg[i], border=edge.lab.border[i], lwd=edge.lab.lwd[i])
			text(edgexy[i,,drop=FALSE], label=edge_labels[i], cex=edge.lab.cex[i], font=edge.lab.font[i], adj=edge.lab.adj[[i]], col=edge.lab.col[i])
		}
	}

# clade lines
	if (isTRUE(show.clade.lines) | isTRUE(show.clade.label)) {
		part <- ape::prop.part(phy)
		labels <- attr(part, "labels")
		nodes <- intersect(phy$edge[backgr,2], phy$edge[!backgr,1])
		units <- vector("list", length(nodes))
		for (i in seq_along(units)) {
			off <- phy$edge[phy$edge[,1] == nodes[i] & !backgr,2]
			units[[i]] <- c(phy$tip.label[off[off <= ntip]], unlist(lapply(part[off[off > ntip] - ntip], function(p) labels[p])))
		}
		names(units) <- sapply(units, function(x) unique(clades[match(x, tips)]))
		clade.lin.col <- color[names(units)]
		ranges <- lapply(units, function(p) range(xy[match(p,tip_labels),2]))
		overhang <- c(-1, 1) * clade.lin.over * ifelse(ycex > 1, tipstrh / ycex, tipstrh)
		ranges <- Map("+", ranges, rep(list(overhang), length(ranges)))
		if (isTRUE(show.clade.lines)) {
			for (i in seq_along(ranges)) {
				axis(side=4, at=ranges[[i]], labels=FALSE, tick=TRUE, lwd=edge.width, lwd.ticks=0, col=clade.lin.col[i], line=clade.lin.pos)
			}
		}
		if (isTRUE(show.clade.label)) {
			if (is.null(clade_labels)) {
				clade_labels <- names(units)
			} else {
				clade_labels <- clade_labels[names(units)]
				if (any(is.na(clade_labels))) {
					clade_labels <- names(units)
					warning("the provided 'clade_labels' are ignored as they do not contain all clade names")
				}
			}
			if (is.null(clade.lab.col)) {
				clade.lab.col <- clade.lin.col
			} else {
				clade.lab.col <- color[clade_labels]
			}
			for (i in seq_along(ranges)) {
				mtext(text=clade_labels[i], side=4, at=mean(ranges[[i]]), line=clade.lin.pos+1.0*onecharw/0.2, cex=clade.lab.cex, font=clade.lab.font, las=clade.lab.las, col=clade.lab.col[i])
			}
		}
	}

# divergence scale
	if (isTRUE(show.scale)) {
		at <- axTicks(side=1, usr=range(xy[,1]))
		step <- diff(at[1:2])
		at <- mean(at) + c(-0.5, 0.5) * step
		if (length(scale.pos) == 2) {
			at <- at + scale.pos[2]
		}
		onecharh <- strheight(s="0", units="inches", cex=1) * scale.cex
		axis(side=1, at=at, labels=FALSE, lwd=edge.width/2, cex.axis=scale.cex, tick=TRUE, lwd.ticks=0, line=scale.pos[1])
		mtext(text=formatC(step, format="f", digits=1), side=1, at=mean(at), cex=scale.cex, line=scale.pos[1]+0.25*onecharh/0.2)
	}

# tidy up
	par(def.par)

# closing pdf device
	if (isTRUE(device == "pdf")) {
		invisible(dev.off())
	}
	
}



drawPhylogram <- function(phy, xy, lwd=1, col=1) {
	lwd <- rep_len(lwd, nrow(phy$edge))
	col <- rep_len(col, nrow(phy$edge))
	for (i in seq(nrow(phy$edge))) {
		edg <- xy[phy$edge[i,],]
		lines(matrix(edg[c(1,2,4,4)],2,2), lwd=lwd[i], col=col[i])
		lines(matrix(edg[c(1,1,3,4)],2,2), lwd=lwd[i], col=col[i])
	}
}

drawCladogram <- function(phy, xy, lwd=1, col=1) {
	lwd <- rep_len(lwd, nrow(phy$edge))
	col <- rep_len(col, nrow(phy$edge))
	for (i in seq(nrow(phy$edge))) {
		lines(xy[phy$edge[i,],], lwd=lwd[i], col=col[i])
	}
}

drawEllipse <- function(x, y, size, w, h, r=0, res=200, col, border, lwd) {
	rs <- seq(0, 2 * pi, len=res)
	pts <- cbind(0.5 * w * cos(rs), 0.5 * h * sin(rs))
	rot <- matrix(c(cos(r), sin(r), -sin(r), cos(r)), 2, 2)
	pts <- size * pts %*% rot + rep(1, res) %*% t(c(x,y))
	polygon(pts, col=col, border=border, lwd=lwd)
}


# Function to dynamically plot text without nested scaling
plot_text <- function(x, y, label, desired_width, font, adj) {
  # Set the text size based on the desired width and reset par("cex") each time
  par(cex = 1)  # Reset cex to default
  cex_adjusted <- desired_width / strwidth(label, cex = 1, font=font)  # Calculate cex with baseline
  # Plot text with adjusted cex
  text(x, y, label, cex = cex_adjusted, font=font, adj=adj)
}

