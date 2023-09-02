# depends on: ape

# inputs:
#	phy - object of class "phylo"
#	rescale - whether to rescale branch lengths
#	bw - bandwidth function
#	kernel - function used in kernel density estimate
#	consp, heterosp - not implemented yet

branchcutting <- function(phy, rescale=TRUE, bw="mdif", kernel="epanechnikov") {

	rescale_branchlengths_rooted <- function(phy) {
		nodes <- lapply(seq(max(phy$edge)), function(i) which(phy$edge[,1] == i))
		deg <- sapply(nodes, function(x) sum(phy$edge.length[x]))
		phy$edge.length <- phy$edge.length * deg[phy$edge[,1]]
		return(phy)
	}

	rescale_branchlengths_unrooted <- function(phy) {
		nodes <- lapply(seq(max(phy$edge)), function(i) which(phy$edge[,1] == i))
		deg <- sapply(nodes, function(x) sum(phy$edge.length[x]))
		phy$edge.length <- phy$edge.length * deg[phy$edge[,1]] * deg[phy$edge[,2]]
		return(phy)
	}

	calculate_loss <- function(phy) {
		edge <- phy$edge
		ntip <- min(edge[,1]) - 1
		noff <- as.numeric(edge[,2] <= ntip)
		last <- which(noff == 1)
		while(any(noff == 0)) {
			term <- last[duplicated(edge[last,1])]
			term <- which(edge[,1] %in% edge[term,1])
			kept <- setdiff(last, term)
			last <- which(edge[,2] %in% edge[term,1])
			new <- tapply(noff[term], edge[term,1], sum)
			noff[last] <- new[order(edge[last,2])]
			last <- c(last, kept)
		}
		npath <- noff * (ntip - noff)
		loss <- npath * phy$edge.length / choose(ntip, 2)
		return(loss)
	}
	
	if (is.character(phy)) {
		if (grepl("^\\(", phy)) {
			phy <- ape::read.tree(text=phy)
		} else {
			phy <- ape::read.tree(phy)
		}
	}
	if (isTRUE(rescale)) {
		rescaling <- list(rescale_branchlengths_unrooted, rescale_branchlengths_rooted)
		rescale_branchlengths <- rescaling[[ape::is.rooted(phy) + 1]]
		rphy <- rescale_branchlengths(phy)
		loss <- calculate_loss(rphy)
	} else {
		rphy <- NULL
		loss <- calculate_loss(phy)
	}
	
	ord <- order(loss)
	loss <- loss[ord]

	if (bw == "mdif") {
		bw <- mean(diff(loss))
	}
		
	dens <- suppressWarnings(density(loss, bw=bw, kernel=kernel, cut=0))
	zero <- sqrt(.Machine$double.eps)
	nonzeros <- which(dens$y <= zero)
	if (length(nonzeros) > 0) {
		bpoint <- dens$x[min(nonzeros)]
	} else {
		bpoint <- max(loss)
	}
	signif <- which(loss > bpoint)

	if (length(signif) > 0) {
		retained <- ord[signif]
		collapsed <- phy
		collapsed$edge.length <- ifelse(seq_along(ord) %in% retained, phy$edge.length, 0)
		retnodes <- unique(sort(phy$edge[retained,]))
		if (ape::is.rooted(phy)) {
			retnodes <- c(retnodes, length(phy$tip.label) + 1)
		}
		retpaths <- vector("list", 0.5 * length(retnodes) * (length(retnodes)-1))
		ij <- 1
		for (i in 1:(length(retnodes)-1)) {
			for (j in (i+1):length(retnodes)) {
				retpaths[[ij]] <- ape::nodepath(phy, retnodes[i], retnodes[j])
				ij <- ij + 1
			}
		}
		retnodes <- unique(unlist(retpaths))
		retedges <- phy$edge[,1] %in% retnodes & phy$edge[,2] %in% retnodes
		collapsed$edge.length[retedges] <- phy$edge.length[retedges]
		otus <- unique(lapply(seq_along(phy$tip.label), function(i, dst) names(which(dst[i,] == 0)), dst=cophenetic(collapsed)))
	} else {
		retained <- numeric(0)
		collapsed <- phy
		otus <- list(phy$tip.label)
	}

	result <- list(phy=phy, otus=otus, notus=length(otus), loss=loss[match(seq_along(loss),ord)], dens=dens, retained=retained, collapsed=collapsed, rescaled=rphy)
	class(result) <- "bcut"

	return(result)
}



plot.bcut <- function(bcut, which=c("tree","otus","loss","dens"), edge.color, edge.width, pt.cex, pt.col=1, ...) {
	if (missing("edge.color")) {
		edge.color <- 1
	}
	if (missing("edge.width")) {
		edge.width <- 1
	}
	if (missing("pt.cex")) {
		pt.cex <- 1
	}
	ord <- order(bcut$loss)
	which <- which[1]
	if (which == "tree") {
		col <- (seq_along(bcut$loss) %in% bcut$retained) + 1
		col[col == 1 & bcut$collapsed$edge.length > 0] <- 3
#		col[col == 2 & bcut$collapsed$edge.length == 0] <- 4
		if (is.rooted(bcut$phy)) {
			plot(bcut$phy, type="phylogram", edge.color=col, edge.width=edge.width, ...)
		} else {
			plot(bcut$phy, type="unrooted", edge.color=col, edge.width=edge.width, ...)
		}
	}
	if (which == "otus") {
		find_mrca <- function(x, phy) ifelse(length(x) == 1, match(x, phy$tip.label), ape::getMRCA(phy, tip=x))
		size <- sapply(bcut$otus, length)
		prop <- size / sum(size)
		retained <- which(bcut$collapsed$edge.length > 0)
		anc <- sapply(bcut$otus, find_mrca, phy=bcut$phy)
		nz <- sort(unique(as.numeric(bcut$phy$edge[retained,])))
		dst <- ape::dist.nodes(bcut$collapsed)
		for (i in which(!anc %in% nz)) {
			anc[i] <- nz[which(dst[anc[i],nz] == 0)]
		}
		type <- c("unrooted","cladogram")[ape::is.rooted(bcut$phy) + 1]

		plot(bcut$phy, type=type, show.tip.label=FALSE, edge.color="transparent", edge.width=edge.width)#, ...)
#		xy <- ape::plotPhyloCoor(bcut$phy, type=type, direction="rightwards", show.tip.label=FALSE, edge.width=edge.width, use.edge.length=TRUE)
		lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
		xy <- cbind(lastPP$xx, lastPP$yy)
		xspl <- vector("list", length(anc))
		ultrametr <- ape::is.ultrametric(bcut$phy, tol=0.0001)
		for (i in seq_along(anc)) {
			tips <- match(bcut$otus[[i]], bcut$phy$tip.label)
			nodes <- lapply(tips, function(x) nodepath(bcut$phy, x, anc[i]))
			nodes <- unique(unlist(nodes))
			nodes <- nodes[chull(xy[nodes,])]
			if (length(nodes) >= 3) {
				xynodes <- scale(xy[nodes,], scale=FALSE)
				xynodes <- 1.1 * xynodes + rep(1, length(nodes)) %*% t(attr(xynodes,"scaled:center"))
				spl <- c(-0.5, 0)[ultrametr + 1]
				xspl[[i]] <- graphics::xspline(x=xynodes[,1], y=xynodes[,2], s=rep(spl, length(nodes)), open=FALSE, draw=FALSE)
			}
		}
		if (ultrametr) {
			xmax <- max(unlist(lapply(xspl, "[[", "x")))
			xtip <- max(xy[ape::Ntip(bcut$phy),1])
			for (i in which(!sapply(xspl, is.null))) {
				xspl[[i]]$x[xspl[[i]]$x >= xtip] <- xmax
			}
		}
		for (i in which(!sapply(xspl, is.null))) {
			x <- data.frame(x=xspl[[i]]$x, y=xspl[[i]]$y)
			ch <- grDevices::chull(x)
			xspl[[i]] <- x[c(ch,ch[1]),]
		}		
		for (i in which(!sapply(xspl, is.null))) {
			polygon(xspl[[i]], col="grey80", border="transparent")
		}
		for (i in seq_along(bcut$phy$edge.length)) {
			lines(xy[bcut$collapsed$edge[i,],], lwd=edge.width, col=edge.color)
		}
		points(xy[anc,,drop=FALSE], pch=16, cex=pt.cex, col=pt.col)
	}
	if (which == "loss") {
		col <- (seq_along(bcut$loss) %in% bcut$retained) + 1
		col[col == 1 & bcut$collapsed$edge.length > 0] <- 3
		col <- col[order(bcut$loss)]
		loss <- sort(bcut$loss)
		plot(seq_along(loss), loss, type="n", xlab="Rank", ylab="Loss", cex.axis=1.2, cex.lab=1.3, ...)
		subs <- col == "green3"
		points(seq_along(loss)[!subs], loss[!subs], pch=16, col=col[!subs], ...)
		points(seq_along(loss)[subs], loss[subs], pch=16, col=col[subs], ...)
	}
	if (which == "dens") {
		dens <- bcut$dens
		zero <- sqrt(.Machine$double.eps)
		bpoint <- dens$x[min(which(dens$y <= zero))]
		signif <- bcut$loss > bpoint
		plot(dens, xlab="Loss", cex.axis=1.2, cex.lab=1.3, ...)
		if (any(signif)) {
			thr <- min(bcut$loss[signif])
			thr <- max(dens$x[dens$x < thr & dens$y <= zero])
			abline(v=thr, col="red")
		}
	}
}


summary.bcut <- function(bcut) {
	cat("Object of class 'bcut':\n")
	cat(paste(paste("\tNo. of OTUs:", length(bcut$otus)), "\n", sep=""))
	cat(paste(paste("\tNo. of tree tips:", length(bcut$phy$tip.label)), "\n", sep=""))
	cat(paste(paste("\tNo. of retained branches:", length(bcut$retained)), "\n", sep=""))
}

							
print.bcut <- function(bcut) {
	k <- length(bcut$otus)
	nn <- sapply(bcut$otus, length)
	o <- ifelse(k == 1, "OTU", "OTUs")

	cat(paste(paste("Object of class 'bcut' with", k, "delimited", o), ":\n", sep=""))
	for (i in seq(k)) {
		six <- head(bcut$otus[[i]])
		if (nn[i] > 6)
			six <- c(six, "...")
		otu <- paste(paste("\t", paste("OTU", i), ":", sep=""), paste(six, collapse=", "))
		cat(paste(otu, "\n", sep=""))
	}
}


write.bcut <- function(bcut, file=NULL, return=TRUE) {
	df <- data.frame(ID=unlist(bcut$otus), OTU=rep(seq_along(bcut$otus), sapply(bcut$otus, length)), stringsAsFactors=FALSE)
	if (!is.null(file)) {
		write.table(df, file, row.names=FALSE, quote=FALSE, sep="\t")
	}
	if (isTRUE(return)) {
		attr(df, "nsp") <- length(bcut$otus)
		return(df)
	}
}

