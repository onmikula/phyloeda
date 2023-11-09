rescale_branchlengths <- function(phy) {
	nodes <- lapply(seq(max(phy$edge)), function(i) which(phy$edge[,1] == i))
	deg <- sapply(nodes, function(x) sum(phy$edge.length[x]))
	if (ape::is.rooted(phy)) {
		phy$edge.length <- phy$edge.length * deg[phy$edge[,1]]		
	} else {
		phy$edge.length <- phy$edge.length * deg[phy$edge[,1]] * deg[phy$edge[,2]]		
	}
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


find_backbone <- function(phy, nodes) {	
	if (ape::is.rooted(phy)) {
		edge <- phy$edge
		last <- edge[,2] %in% nodes
		backbone <- last
		while(any(last)) {
			backbone <- backbone | last
			last <- edge[,2] %in% edge[last,1]
		}
		backbone <- backbone | edge[,1] %in% edge[backbone,1]  & !backbone
	} else {
		edge <- apply(apply(phy$edge, 1, sort), 2, paste, collapse="-")
		backbone <- rep(FALSE, length(edge))
		terminal <- sort(nodes)
		while (length(terminal) > 1) {
			path <- ape::nodepath(phy, terminal[1], terminal[length(terminal)])
			terminal <- setdiff(terminal, path)
			path <- cbind(path[-1], path[-length(path)])
			path <- apply(apply(path, 1, sort), 2, paste, collapse="-")
			backbone <- backbone | edge %in% path
		}
		if (length(terminal) == 1) {
			path <- ape::nodepath(phy, terminal, setdiff(nodes, terminal)[1])
			path <- cbind(path[-1], path[-length(path)])
			path <- apply(apply(path, 1, sort), 2, paste, collapse="-")
			backbone <- backbone | edge %in% path			
		}
		bbnodes <- unique(as.numeric(phy$edge[backbone,]))
		backbone <- phy$edge[,1] %in% bbnodes | phy$edge[,2] %in% bbnodes
	}
	return(backbone)
}


find_terminals <- function(phy, nodes, tips) {
	if (mode(tips) != "numeric") {
		tips[,] <- match(tips, phy$tip.label)
		mode(tips) <- "numeric"
	}
	edge <- phy$edge
	nedge <- nrow(edge)
	ntip <- length(phy$tip.label)
	terminals <- rep(FALSE, nedge)
	for (i in seq_along(nodes)) {
		path1 <- ape::nodepath(phy, nodes[i], tips[i,1])
		path2 <- ape::nodepath(phy, nodes[i], tips[i,2])
		paths <- rbind(cbind(path1[-length(path1)], path1[-1]), cbind(path2[-length(path2)], path2[-1]))
		terminals <- terminals | duplicated(rbind(edge, paths), fromLast=TRUE)[1:nedge]
		offshoots <- !terminals & edge[,1] %in% paths
		while (any(offshoots)) {
			terminals <- terminals | offshoots
			offshoots <- !terminals & edge[,1] %in% edge[offshoots,2]
		}
	}
	return(terminals)
}


find_outliers <- function(loss, bw, kernel, backbone) {
	dens <- suppressWarnings(density(loss[!backbone], bw=bw, kernel=kernel, cut=0))
	zero <- sqrt(.Machine$double.eps)
	nonzeros <- which(dens$y <= zero)
	if (length(nonzeros) > 0) {
		bpoint <- dens$x[min(nonzeros)]
	} else {
		bpoint <- max(loss)
	}
	outliers <- loss > bpoint
	result <- list(outliers=outliers, dens=dens, bpoint=bpoint)
	return(result)	
}


# depends on: ape

# inputs:
#	phy - object of class "phylo"
#	heterosp, consp - 
#	rescale - whether to rescale branch lengths
#	bw - bandwidth function
#	kernel - function used in the kernel density estimation
#	full.output - whether the full output is returned

bcut_core <- function(phy, heterosp=NULL, consp=NULL, rescale=TRUE, bw="mdif", kernel="epanechnikov", full.output=TRUE) {
	if (is.character(phy)) {
		if (grepl("^\\(", phy)) {
			phy <- ape::read.tree(text=phy)
		} else {
			phy <- ape::read.tree(phy)
		}
	}

	if (isTRUE(rescale)) {
		rphy <- rescale_branchlengths(phy)
		loss <- calculate_loss(rphy)
	} else {
		rphy <- NULL
		loss <- calculate_loss(phy)
	}
	
	ntip <- length(phy$tip.label)
	nedge <- nrow(phy$edge)
	if (!is.null(heterosp)) {
		if (length(heterosp) == 1) {
			heterosp <- read.delim(heterosp)
		}
		maxdeepanc <- sapply(seq(nrow(heterosp)), function(i) ape::getMRCA(phy, tip=heterosp[i,]))
		hetconstrain <- find_backbone(phy, nodes=maxdeepanc)
	} else {
		hetconstrain <- rep(FALSE, nedge)
	}
	if (!is.null(consp)) {
		if (length(consp) == 1) {
			consp <- read.delim(consp)
		}
		mindeepanc <- sapply(seq(nrow(consp)), function(i) ape::getMRCA(phy, tip=consp[i,]))
		conconstrain <- find_terminals(phy, nodes=mindeepanc, tips=consp)
	} else {
		conconstrain <- rep(FALSE, nedge)
	}
	if (any(hetconstrain & conconstrain)) {
		stop("Heterospecific and conspecific constraints are in conflict.")
	}
	
	if (bw == "mdif") {
		bw <- mean(diff(sort(loss)))
	}
	outliers <- find_outliers(loss, bw=bw, kernel=kernel, backbone=hetconstrain)
	dens <- outliers$dens
	bpoint <- outliers$bpoint
	outliers <- outliers$outliers
	retained <- sort(which((outliers | hetconstrain) & !conconstrain))
	if (length(retained) > 0) {
		retnodes <- unique(sort(phy$edge[retained,]))
		backbone <- find_backbone(phy, nodes=retnodes)
		ancestors <- setdiff(phy$edge[backbone, 2], phy$edge[backbone, 1])
		singletons <- ancestors <= ntip
		otus <- vector("list", length(ancestors))		
		otus[singletons] <- phy$tip.label[ancestors[singletons]]
		otus[!singletons] <- lapply(ancestors[!singletons], function(i) ape::extract.clade(phy, node=i)$tip.label)
	} else {
		otus <- list(phy$tip.label)
	}

	if (isTRUE(full.output)) {
		otus <- data.frame(ID=unlist(otus), OTU=rep(seq_along(otus), sapply(otus, length)), stringsAsFactors=FALSE)
		otus$OTU <- paste0("BCUT", formatC(otus$OTU, format="d", flag=0, width=nchar(max(otus$OTU))))
		criteria <- data.frame(backbone=backbone, outliers=outliers, heterosp=hetconstrain, consp=conconstrain)
		output <- list(otus=otus, phy=phy, loss=loss, dens=dens, bpoint=bpoint, backbone=which(backbone), criteria=criteria, notus=length(unique(otus[,2])), rescaled=rphy)
		class(output) <- "bcut"
	} else {
		otus <- lapply(otus, sort)
		output <- otus
	}
	return(output)
}


branchcutting <- function(phy, psample, heterosp=NULL, consp=NULL, rescale=TRUE, bw="mdif", kernel="epanechnikov") {

	if (!missing(phy)) {
		if (is.character(phy) & length(phy) == 1) {
			phy <- ape::read.tree(phy)
		}
	} else {
		phy <- NULL
	}
	if (is.list(phy) & !inherits(phy, "phylo") & missing(psample)) {
		psample <- phy
		phy <- NULL
	} else if (is.list(phy) & !inherits(phy, "phylo") & !missing(psample)) {
		phy <- NULL
		warning("cannot work with 'phy' as it is not in 'phylo' format")
	} else if (missing(psample)) {
		psample <- NULL
	}

	if (!is.null(phy)) {
		bcut <- bcut_core(phy=phy, heterosp=heterosp, consp=consp, rescale=rescale, bw=bw, kernel=kernel, full.output=TRUE)		
	}
	
	if (!is.null(psample)) {
		ps <- lapply(psample, bcut_core, heterosp=heterosp, consp=consp, rescale=rescale, bw=bw, kernel=kernel, full.output=FALSE)
		ps <- lapply(unlist(ps, recursive=FALSE), sort)
		pp <- table(sapply(ps, paste, collapse="-")) / length(psample)
		pp <- pp[order(pp, names(pp), method="radix", decreasing=c(TRUE, FALSE))]
		otus <- strsplit(names(pp), "-")
		names(otus) <- paste0("OTU", formatC(seq_along(otus), format="d", flag=0, width=nchar(length(otus))))
		tips <- psample[[1]]$tip.label
		mat <- matrix(0, length(otus), length(otus))
		for (i in 1:(length(otus) - 1)) {
			for (j in (i+1):length(otus)) {
				mat[i,j] <- mat[i,j] <- as.numeric(length(intersect(otus[[i]], otus[[j]])) == 0)
			}
		}
		mcdelim <- 1
		complete <- length(setdiff(tips, unlist(otus[mcdelim]))) == 0
		while (complete == FALSE) {		
			mcdelim <- c(mcdelim, min(which(colSums(mat[mcdelim,,drop=FALSE] == 1) == length(mcdelim))))
			complete <- length(setdiff(tips, unlist(otus[mcdelim]))) == 0
		}
		otudf <- data.frame(OTU=names(otus), PP=as.numeric(pp), MC=seq_along(otus) %in% mcdelim, row.names=NULL)
		map <- matrix(0, length(otus), length(tips), dimnames=list(NULL, tips))
		for (i in seq_along(otus)) {
			map[i,otus[[i]]] <- 1
		}
		psample <- list(pp=otudf, classification=map)
		if (!is.null(phy)) {
			otus_bcut <- lapply(split(bcut$otus[,1], bcut$otus[,2]), sort)
			bcut$pp <- data.frame(OTU=names(otus_bcut), PP=as.numeric(pp[match(otus_bcut, otus)]))
			bcut$psample <- psample
		}
	}
	
	if (!is.null(phy)) {
		output <- bcut
		class(output) <- "bcut"
	} else {
		mc_otus <- data.frame(ID=unlist(otus[otudf$MC]), OTU=rep(names(otus[otudf$MC]), sapply(otus[otudf$MC], length)), stringsAsFactors=FALSE, row.names=NULL)
		output <- list(otus=mc_otus, psample=psample)
	}
	return(output)
	
}



summary.bcut <- function(bcut) {
	cat("Object of class 'bcut':\n")
	cat(paste(paste("\tNo. of OTUs:", bcut$notus), "\n", sep=""))
	cat(paste(paste("\tNo. of tree tips:", length(bcut$phy$tip.label)), "\n", sep=""))
	cat(paste(paste("\tNo. of retained branches:", length(bcut$backbone)), "\n", sep=""))
}

							
print.bcut <- function(bcut) {
	k <- bcut$notus
	nn <- table(bcut$otus[,2])
	o <- ifelse(k == 1, "OTU", "OTUs")

	cat(paste(paste("Object of class 'bcut' with", k, "delimited", o), ":\n", sep=""))
	for (i in seq(k)) {
		six <- head(bcut$otus[,1][bcut$otus[,2] == names(nn)[i]])
		if (nn[i] > 6)
			six <- c(six, "...")
		otu <- paste(paste("\t", paste("OTU", i), ":", sep=""), paste(six, collapse=", "))
		cat(paste(otu, "\n", sep=""))
	}
}


write.bcut <- function(bcut, file=NULL, return=TRUE, col.names=c("ID", "OTU"), taxonomy=NULL) {
	otus <- bcut$otus
	if (!is.null(taxonomy)) {
		taxonomy <- taxonomy[match(otus[,1], taxonomy[,1]), 2]
		tab <- table(otus[,2], taxonomy)
		rows <- rowSums(tab != 0) == 1
		cols <- colSums(tab != 0) == 1
		rows <- names(which(rowSums(tab[rows,cols,drop=FALSE] != 0) == 1))
		cols <- names(which(colSums(tab[rows,cols,drop=FALSE] != 0) == 1))
		tab <- tab[rows,cols,drop=FALSE]
		map <- setNames(colnames(tab)[apply(tab > 0, 1, which)], rownames(tab))
		named <- otus[,2] %in% names(map)
		otus[named,2] <- map[otus[named,2]]
	}
	otus <- setNames(otus, col.names)
	if (!is.null(file)) {
		write.table(df, file, row.names=FALSE, quote=FALSE, sep="\t")
	}
	if (isTRUE(return)) {
		attr(df, "notus") <- bcut$notus
		return(df)
	}
}



plot.bcut <- function(bcut, which=c("tree","otus","loss","dens"), edge.width=2, show.tip.label=FALSE, tip.cex=0.7, pt.cex=2, show.legend=FALSE, ...) {

	draw_phylogram <- function(phy, xy, lwd, col, ...) {
		lwd <- rep_len(lwd, nrow(phy$edge))
		col <- rep_len(col, nrow(phy$edge))
		for (i in seq(nrow(phy$edge))) {
			edg <- xy[phy$edge[i,],]
			lines(matrix(edg[c(1,2,4,4)],2,2), lwd=lwd[i], col=col[i], ...)
			lines(matrix(edg[c(1,1,3,4)],2,2), lwd=lwd[i], col=col[i], ...)
		}
	}
	draw_unrooted <- function(phy, xy, lwd, col, ...) {
		lwd <- rep_len(lwd, nrow(phy$edge))
		col <- rep_len(col, nrow(phy$edge))
		for (i in seq(nrow(phy$edge))) {
			edg <- xy[phy$edge[i,],]
		}
	}

	which <- which[1]

	if (which == "tree") {
		col <- ifelse(bcut$criteria$backbone, 3, 1)
		col <- ifelse(bcut$criteria$outliers, 2, col)
		col <- ifelse(!bcut$criteria$outliers & bcut$criteria$heterosp, 4, col)
		col <- ifelse(bcut$criteria$outliers & bcut$criteria$consp, 8, col)
		type <- ifelse(ape::is.rooted(bcut$phy), "phylogram", "unrooted")
		draw_tree <- list(phylogram=draw_phylogram, unrooted=draw_unrooted)[[type]]
		xy <- ape::plotPhyloCoor(bcut$phy, type=type, direction="rightwards", use.edge.length=TRUE)
		xlim <- range(xy[,1])
		if (.Device == "null device") {
			device <- match.fun(ifelse(.Platform$OS.type == "unix", "quartz", "x11"))
			device(width=7, height=7)
			par(mai=c(0.22, 0.22, 0.22, 0.22))
		}
		if (isTRUE(show.tip.label)) {
			ntip <- length(phy$tip.label)
			offset <- 0.25 * mean(phy$edge.length[phy$edge[,2] <= ntip])
			xlim[2] <- xlim[2] + offset
			strw <- tip.cex * max(strwidth(phy$tip.label, units="inches", font=3))
			strw <- strw * 1.08 * diff(xlim) / (par("fin")[2] - sum(par("mai")[c(2,4)]))
			xlim[2] <- xlim[2] + strw
		}
		plot(xy, type="n", bty="n", axes=FALSE, ann=FALSE, xlim=xlim)
		draw_tree(phy, xy, lwd=edge.width, col=col)
		if (isTRUE(show.tip.label)) {
			xytips <- xy[seq(ntip),]
			xytips[,1] <- xytips[,1] + offset
			text(xytips, labels=phy$tip.label, cex=tip.cex, font=3, adj=c(0,0.25))
		}
	}

	if (which == "otus") {
		find_mrca <- function(x, phy) ifelse(length(x) == 1, phy$edge[match(match(x,phy$tip.label),phy$edge[,2]),1], ape::getMRCA(phy, tip=x))
		draw_ellipse <- function(size, w, h) {
			rs <- seq(0, 2 * pi, len=200)
			pts <- scale(cbind(0.5 * w * cos(rs), 0.5 * h * sin(rs)), scale=FALSE)
			return(pts %*% diag(rep(size / w, 2)))
		}
		backbone <- bcut$backbone
		anc <- sapply(split(bcut$otus[,1], bcut$otus[,2]), find_mrca, phy=bcut$phy)
		type <- c("unrooted","cladogram")[ape::is.rooted(bcut$phy) + 1]
		xy <- ape::plotPhyloCoor(bcut$phy, type=type, direction="rightwards", edge.width=edge.width, use.edge.length=TRUE)
		offset <- mean(bcut$phy$edge.length[bcut$phy$edge[,2] <= length(bcut$phy$tip.label)])
		circle <- draw_ellipse(size=offset, w=diff(range(xy[,1])), h=diff(range(xy[,2])))
		hulls <- vector("list", length(anc))
		for (i in seq_along(anc)) {
			tips <- match(bcut$otus[bcut$otus[,2] == names(anc)[i],1], bcut$phy$tip.label)
			nodes <- unique(unlist(lapply(tips, function(x) ape::nodepath(bcut$phy, x, anc[i]))))
			nodes <- nodes[grDevices::chull(xy[nodes,])]
			hulls[[i]] <- lapply(nodes, function(n, circle) circle + rep(1, 200) %*% t(xy[n,]), circle=circle)
			hulls[[i]] <- do.call(rbind, hulls[[i]])
			hulls[[i]] <- hulls[[i]][grDevices::chull(hulls[[i]]),]
		}
		xylim <- do.call(rbind, lapply(c(list(xy), hulls), function(x) apply(x, 2, range)))
		xylim <- apply(xylim, 2, range)
		if (.Device == "null device") {
			device <- match.fun(ifelse(.Platform$OS.type == "unix", "quartz", "x11"))
			device(width=7, height=7)
			par(mai=c(0.22, 0.22, 0.22, 0.22))
		}
		plot(xy, type="n", bty="n", axes=FALSE, ann=FALSE, xlim=xylim[,1], ylim=xylim[,2])
		for (i in seq_along(anc)) {
			polygon(hulls[[i]], col="grey80", border="transparent")
		}
		for (i in seq(nrow(bcut$phy$edge))) {
			lines(xy[bcut$phy$edge[i,],], lwd=edge.width, col="black")
		}
		single <- table(bcut$otus[,2])[names(anc)] == 1
		if (any(single)) {
			off <- bcut$otus[,2] %in% names(anc)[single]
			off <- match(bcut$otus[off,1], bcut$phy$tip.label)
			points(xy[off,,drop=FALSE], pch=16, cex=pt.cex, col="black")
		}
		if (any(!single)) {
			points(xy[anc[!single],,drop=FALSE], pch=16, cex=pt.cex, col="black")
		}
	}
	
	if (which == "loss") {
		col <- ifelse(bcut$criteria$backbone, 3, 1)
		col <- ifelse(bcut$criteria$outliers, 2, col)
		col <- ifelse(!bcut$criteria$outliers & bcut$criteria$heterosp, 4, col)
		col <- ifelse(bcut$criteria$outliers & bcut$criteria$consp, 8, col)
		col <- col[order(bcut$loss)]
		loss <- sort(bcut$loss)
		if (.Device == "null device") {
			device <- match.fun(ifelse(.Platform$OS.type == "unix", "quartz", "x11"))
			device(width=7, height=7)
		}
		par(mai=c(1.02, 0.92, 0.42, 0.42))
		plot(seq_along(loss), loss, type="n", xlab="Rank", ylab="Loss", main="", cex.axis=1.2, cex.lab=1.3, ...)
		points(which(col == 1), loss[col == 1], pch=16, col=1, ...)
		points(which(col == 3), loss[col == 3], pch=16, col=3, ...)
		points(which(col == 4), loss[col == 4], pch=16, col=4, ...)
		points(which(col == 8), loss[col == 8], pch=16, col=8, ...)
		points(which(col == 2), loss[col == 2], pch=16, col=2, ...)
	}
	
	if (which == "dens") {
		if (.Device == "null device") {
			device <- match.fun(ifelse(.Platform$OS.type == "unix", "quartz", "x11"))
			device(width=7, height=7)
		}
		par(mai=c(1.02, 0.92, 0.42, 0.42))
		plot(bcut$dens, xlab="Loss", main="", cex.axis=1.2, cex.lab=1.3, lwd=2, ...)
		if (any(bcut$criteria$outliers)) {
			abline(v=bcut$bpoint, col="red", lwd=2, lty=3, ...)
		}
	}
}



