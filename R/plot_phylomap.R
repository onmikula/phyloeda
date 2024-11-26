#' Phylogeographic map.
#' 
#' @description
#' Creates a simple phylogeographic map.
#' 
#' @param info a data frame with individual identifier, classification into phylogeographic lineages
#'   and geographical coordinates, whose names are partial matches of "latitude", "longitude" or "x", "y").
#' @param ind,lin optional character strings with column names including individual IDs (e.g. "ID")
#'   and classification to phylogeographic lineage (e.g. "Species"), respectively. If not supplied,
#'   these two columns are supposed to be the first and the second one, respectively.
#' @param gps optional data frame with geographical coordinates. It must contain the same individual identifier 
#'   as used in `info` and geographical coordinates named as explained for `info`. If `gps` is supplied, it supersedes
#'   the information in `info`. If only a subset of IDs is matching those in `info`, only this subset is plotted.
#' @param subset logical, numeric or character; indicates subset of rows in info to be plotted. If logical, it must be
#'   of the same length as it is the number of rows in info, if character, it contains individual IDs.
#' @param color a named vector of point colors, names must correspond to lineage labels.
#'   If missing, a default palette is used.
#' @param border color of point borders. If `border == "auto"`, the borders are identical to point fillings.
#' @param size size of points specified as a percentage of the lesser of latitude / longitude ranges.
#' @param map map layer serving as background instead of country borders from `maps` package.
#' @param device either `"quartz"`, `"x11"` or `"pdf"` to indicate the type of graphical device.
#'   If `NULL`, the objects are plotted into the current device.
#' @param file a name of pdf file if `device == "pdf"`.
#' @param width width of graphical device in inches.
#' @param height height of graphical device in inches.
#' @param mai size of outer margins in inches, recycled if necessary.
#' @param xlim,ylim Limits of x and y axes.
#' @param axes logical, whether to display axes.
#' @param cex.axis magnification of axis annotation.
#' @param legend logical, whether to show the legend (on a margin, with in-built parameters).
#' @param expr an expression (or a list of them) allowing to add other elements (e.g. customized legend).
#' @export

plot_phylomap <- function(info, ind, lin, gps, subset, color, border="black", size=4, lwd=1, map=NULL, device, width=7, height=7, file=NULL, mai=0.02, xlim, ylim, axes=TRUE, cex.axis=1.25, legend=TRUE, expr) {
		
	info <- as.data.frame(info)
	if (missing(ind)) {
		ind <- names(info)[1]
	}
	if (missing(lin)) {
		lin <- names(info)[2]
	}
	xy <- pmatch(c("lon","lat"), tolower(names(info)))
	if (any(is.na(xy))) {
		if (missing(gps)) {
			stop("geographical coordinates are missing or their names are invalid in 'info' and 'gps' is not supplied")
		}
		xy <- NULL
	} else {
		xy <- names(info)[xy]
	}
	info <- info[,c(ind, lin, xy)]
	if (!missing(subset)) {
		if (is.character(subset)) {
			subset <- info[,1] %in% subset
		} else if (is.numeric(subset)) {
			subset <- seq(nrow(info)) %in% subset
		}
		info <- info[subset,]
	}
	if (is.null(xy)) {
		gps <- as.data.frame(gps)
		xy <- pmatch(c("lon","lat"), tolower(names(gps)))
		if (any(is.na(xy))) {
			stop("geographical coordinates are missing or their names are invalid in 'gps'")
		} else {
			xy <- names(gps)[xy]
		}
		if (!any(grepl(ind, names(gps)))) {
			if (!any(grepl(ind, names(gps), ignore.case=TRUE))) {
				stop("individual identifier is missing or its name is invalid in 'gps'")
			} else {
				gps <- setNames(gps[,c(grep(ind, names(gps), ignore.case=TRUE, value=TRUE), xy)], c(ind, xy))
			}
		} else {
			gps <- gps[,c(ind, xy)]
		}
		info <- info[info[,ind] %in% gps[,ind],]
		info <- data.frame(info, gps[match(info[,ind], gps[,ind]),xy])
	}
	
	xy <- info[,xy]
	sp <- info[,lin]
	L <- sort(unique(sp))
	K <- length(L)
	if (missing(color)) {
		palette <- c(Red="#e6194b", Green="#3cb44b", Yellow="#ffe119", Blue="#0082c8", Orange="#f58231", Purple="#911eb4", Cyan="#46f0f0", Magenta="#f032e6", Lime="#d2f53c", Pink="#fabebe", Teal="#008080", Lavender="#e6beff", Brown="#aa6e28", Beige="#fffac8", Maroon="#800000", Mint="#aaffc3", Olive="#808000", Coral="#ffd8b1", Navy="#000080", Grey="#808080")
		color <- setNames(palette[seq(K)], L)
	}

	if (missing(xlim)) {
		xlim <- range(xy[,1], na.rm=TRUE)
	}
	if (missing(ylim)) {
		ylim <- range(xy[,2], na.rm=TRUE)
	}
	circ <-  0.04 * min(diff(xlim), diff(ylim)) * 0.5 * cos(0)
	xlim <- xlim + c(-1, 1) * circ
	ylim <- ylim + c(-1, 1) * circ
	loc <- unique(xy)
	sites <- vector("list", nrow(loc))
	for (i in seq_along(sites)) {
		sites[[i]] <- which(xy[,1] == loc[i,1] & xy[,2] == loc[i,2])
	}
	prop <- lapply(sites, function(x, sp) table(sp[x]), sp=sp)
	ord <- order(sapply(prop, length), sapply(prop, sum))
	loc <- loc[ord,]
	prop <- prop[ord]
	
	# plotting
	if (missing(device)) {
		device <- ifelse(.Platform$OS.type == "unix", "quartz", "x11")
	}
	if (isTRUE(device == "pdf") & missing(file)) {
		file <- "phylo_map.pdf"
	}
	if (isTRUE(device == "pdf")) {
		pdf(file, width=width, height=height)
	} else if (!is.null(device)) {
		match.fun(tolower(device))(width=width, height=height)	
	}
#	opar <- par(no.readonly = TRUE)
	mai <- rep_len(mai, 4)
	if (isTRUE(legend)) {
		mai[4] <- mai[4] + 1.1
		width <- width + 1.1
	}
	if (isTRUE(axes)) {
		onechar <- strheight(s="L", units="inches", cex=1) * cex.axis
		mai[1] <- max(c(mai[1], 3 * onechar + 0.5 * onechar / 0.2))
		mai[2] <- max(c(mai[2], 3 * onechar + 0.5 * onechar / 0.2))
	}
	
	par(mai=mai)
	plot(xlim, ylim, asp=1, type="n", bty="n", axes=FALSE, ann=FALSE)
	if (!is.null(map)) {
		plot(map, asp=1, xlim=xlim, ylim=ylim, add=TRUE)
	} else {
		maps::map(xlim=xlim, ylim=ylim, add=TRUE)
	}

	size <- 0.01 * size * min(diff(xlim), diff(ylim))
	if (border == "auto") {
		border <- lapply(prop, function(x) color[names(x)])
	} else {
		border <- rep(border, length(prop))
	}
	for (i in seq_along(sites)) {
		partCircle(x=loc[i,1], y=loc[i,2], prop=prop[[i]], size=size, res=200, col=color[names(prop[[i]])], border=border[[i]], lwd=lwd)
	}

	if (isTRUE(axes)) {
		line <- 0.5 * onechar / 0.2
		axis(side=1, cex.axis=cex.axis, line=line)
		axis(side=2, cex.axis=cex.axis, line=line)
		title(xlab="Longitude", ylab="Latitude", cex.lab=cex.axis)
	}
	if (isTRUE(legend)) {
		inset <- c(-1 * 1.1 / (width - mai[2] - mai[4]), mai[3] / (width - mai[1] - mai[3]))
		legend(x = "topright", inset=inset, legend=L, pch=21, pt.bg=color[L], pt.lwd=1, xpd=TRUE, bg="white")
	}

	if (!missing(expr)) {
		if (is.expression(expr)) {
			eval(expr)
		} else if (all(sapply(expr, is.expression))) {
			for (i in seq_along(expr)) {
				eval(expr[[i]])
			}
		}		
	}

#	on.exit(par(opar))
# closing pdf device	
	if (isTRUE(device == "pdf")) {
		invisible(dev.off())
	}

}




partCircle <- function(x, y, prop, size, res, col, border, lwd) {
	rs <- seq(0, 2 * pi, len=res)
	pts <- cbind(0.5 * cos(rs), 0.5 * sin(rs))
	pts <- size * pts + rep(1, res) %*% t(c(x,y))
	prop <- prop / sum(prop)
	parts <- round(quantile(seq(nrow(pts)), probs=c(0,cumsum(prop))))
	if (is.null(names(border)) | !all(names(border) %in% names(col))) {
		border <- setNames(rep_len(border, length(col)), names(col))
	}
	if (all(names(prop) %in% names(col)) & !is.null(names(prop))) {
		col <- col[names(prop)]
		border <- border[names(prop)]
	} else {
		col <- rep_len(col, length(prop))
		border <- rep_len(border, length(prop))
	}
	if (length(prop) == 1) {
		polygon(pts, border=border, col=col, lwd=lwd)
	} else {
		for (i in seq_along(prop)) {
			ipts <- rbind(c(x,y), pts[parts[i]:parts[i+1],], c(x,y))
			polygon(ipts, border="transparent", col=col[i], lwd=1)
			lines(ipts[-c(1,nrow(ipts)),], col=border[i], lwd=lwd)
		}	
	}
}

