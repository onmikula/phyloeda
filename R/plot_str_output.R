#' Plotting STRUCTURE output file.
#' 
#' @description
#' Creates barplot of STRUCTURE output.
#'
#' @param out data frame produced by [read_str_output].
#' @param palette character, colors assigned to the clusters in their consecutive order.
#' @param interval logical, whether to display width of probability interval by color transparency.
#' @param ord character or logical, order of individuals in the plot (from left to right).
#'   If `TRUE` (default) it is determined internally, if `FALSE` it is ignored.
#' @param device either `"quartz"`, `"x11"` (= `"X11"`) or `"pdf"`.
#'   If `NULL`, the objects are plotted into the current device.
#' @param file character strigm the name of pdf file if `device == "pdf"`.
#' @param width width of the device in inches.
#' @param height height of the device in inches.
#' @param mai numeric vector, sizes of outer margins (bottom, left, top, right) in inches.
#' @param show.ind.label logical, whether to display individual labels.
#' @param side.ind.label numeric, where to display individual labels (default is `1`, i.e., below the plot).
#' @param cex.axis size of tick labels.
#' @param cex.lab size of tick labels.
#' @param expr expression allowing to add other elements (e.g. legend).
#' @export

plot_str_output <- function(out, palette, interval=FALSE, ord=TRUE, device, file, width=10, height=5, mai=NULL, show.ind.label=FALSE, side.ind.label=1, ylab="Ancestry proportions", main=NULL, line.axis=0.25, cex.axis=1, cex.lab=1, cex.main=1, font.lab=1, expr, ...) {
	N <- nrow(out)
	if (isTRUE(interval) & !any(grepl("low|upp", colnames(out)))) {
		interval <- FALSE
	}
	if (isTRUE(interval)) {
		est <- out[,!grepl("low|upp", colnames(out))]
		low <- out[,grepl("low", colnames(out))]
		upp <- out[,grepl("upp", colnames(out))]
		alpha <- setNames(1 - (upp - low), colnames(est))
		K <- ncol(out) / 3
	} else {
		est <- out[,!grepl("low|upp", colnames(out))]
		alpha <- matrix(1, nrow(est), ncol(est), dimnames=list(rownames(est), colnames(est)))
		K <- ncol(est)
	}
	rankorder <- identical(seq(nrow(out)), sort(unname(ord)))
	nameorder <- identical(sort(rownames(out)), sort(unname(ord)))
	if (!rankorder & !nameorder) {
		if (isTRUE(ord)) {
			ord <- stats::hclust(stats::dist(est), method="average")$order
		} else if (isFALSE(ord)) {
			ord <- seq(N)
		} else {
			ord <- seq(N)
			warning("'ord' argument was set to FALSE")
		}	
	} else if (nameorder) {
		ord <- match(ord, rownames(out))
	}
	est <- est[ord,,drop=FALSE]
	alpha <- alpha[ord,,drop=FALSE]
	csum <- cbind(0, t(apply(est, 1, cumsum)))
	csum <- diag(1/csum[,K+1]) %*% csum

	if (missing(palette)) {
		palette <- c(Red="#e6194b", Green="#3cb44b", Yellow="#ffe119", Blue="#0082c8", Orange="#f58231", Purple="#911eb4", Cyan="#46f0f0", Magenta="#f032e6", Lime="#d2f53c", Pink="#fabebe", Teal="#008080", Lavender="#e6beff", Brown="#aa6e28", Beige="#fffac8", Maroon="#800000", Mint="#aaffc3", Olive="#808000", Coral="#ffd8b1", Navy="#000080", Grey="#808080", White="#FFFFFF", Black="#000000")
	}
	color <- setNames(palette[seq(K)], colnames(est))
	color <- do.call(cbind, lapply(seq(K), function(i) opacity(color[i], alpha[,i])))
	if (is.null(mai)) {
		mai <- c(0.42, 0.82, 0.42, 0.42)
		if (isTRUE(show.ind.label)) {
			mai[c(1,3)] <- c(0.92, 0.32)
		}
	}

	if (missing(device)) {
		if (.Platform$OS.type == "unix") {
			device <- "quartz"
		} else {
			device <- "x11"
		}
	}
	if (isTRUE(device == "pdf")) {
		file <- ifelse(missing(file), "str_barplot.pdf", file)
		pdf(file, width=width, height=height)
	} else if (!is.null(device)) {
		match.fun(tolower(device))(width=width, height=height)	
	}
	par(mai=mai)
	plot(0, 0, xlim=c(0, nrow(out)), ylim=c(0,1), xaxs="i", yaxs="i", xaxt="n", xlab="", ylab=ylab, main=main, type="n", cex.lab=cex.lab, font.lab=font.lab, cex.main=cex.main)
	for (i in seq(N)) {
		x <- rep(i, 5) - c(1,0,0,1,1)
		for (j in seq(K)) {
			xy <- cbind(x, csum[i,j+c(0,0,1,1,0)])
			polygon(xy, border=color[i,j], col=color[i,j])
		}
	}
	if (isTRUE(show.ind.label)) {
		mtext(text=rownames(out)[ord], side=side.ind.label, line=line.axis, at=seq(0.5, nrow(out)-0.5, by=1), cex=cex.axis, las=2, ...)
	}
	if (!missing(expr)) {
		if (is.expression(expr)) {
			eval(expr)
		}
	}
	if (isTRUE(device == "pdf")) {
		invisible(dev.off())
	}
}



#' Custom STRUCTURE palette.
#' 
#' @description
#' Prepares a custom palette for [plot_str_output].
#'
#' @param out data frame produced by [read_str_output].
#' @param info a data frame listing individual labels (1st column) and population labels (2nd column).
#' @param palette a named vector of colors, the names are labels of populations from `info`.
#' @export

custom_str_palette <- function(out, info, palette) {
	trubetskoy <- c(Red="#e6194b", Green="#3cb44b", Yellow="#ffe119", Blue="#0082c8", Orange="#f58231", Purple="#911eb4", Cyan="#46f0f0", Magenta="#f032e6", Lime="#d2f53c", Pink="#fabebe", Teal="#008080", Lavender="#e6beff", Brown="#aa6e28", Beige="#fffac8", Maroon="#800000", Mint="#aaffc3", Olive="#808000", Coral="#ffd8b1", Navy="#000080", Grey="#808080", White="#FFFFFF", Black="#000000")
	K <- attr(out ,"K")
	cluster <- apply(out [,seq(K)], 1, which.max)
	population <- info[match(names(cluster), info[,1]), 2]
	matching <- sort(apply(table(population, cluster), 1, which.max))
	if (any(duplicated(matching))) {
		pal <- trubetskoy[1:K]
		warning("no unique matching, default palette is returned")
	} else if (length(matching) < K) {
		ord <- order(c(matching, setdiff(seq(K), matching)))
		pal <- unique(c(palette[names(matching)], trubetskoy))[ord]
	} else {
		pal <- palette[names(matching)]
	}
	return(pal)
}
