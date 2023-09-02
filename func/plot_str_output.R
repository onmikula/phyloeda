### plot_str_output ###
# out: data frame produced by read_str_output
# palette: colors assigned to the clusters
# interval: whether to display width of probability interval by color transparency
# ord: order of individuals in the plot, if TRUE (default) it is determined internally, if FALSE it is ignored 
# device: "quartz", "x11" (or "X11") or "pdf", if NULL, the objects are plotted into the current device
# file: name of pdf file if device == "pdf"
# width: width of the device in inches
# height: height of the device in inches
# mai: 
# show.ind.label: whether to display individual labels below the plot
# cex.axis: size of tick labels
# cex.lab: size of tick labels
# expr: expression allowing to add other elements (e.g. legend)

plot_str_output <- function(out, palette, interval=FALSE, ord=TRUE, device, file, width=10, height=5, mai=NULL, show.ind.label=FALSE, ylab="Ancestry proportions", cex.axis=1, cex.lab=1, font.lab=1, expr) {
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
	if (!identical(sort(rownames(out)), sort(unname(ord)))) {
		if (isTRUE(ord)) {
			ord <- stats::hclust(stats::dist(est), method="average")$order
		} else if (isFALSE(ord)) {
			ord <- seq(N)
		} else {
			ord <- seq(N)
			warning("'ord' argument was set to FALSE")
		}	
	} else {
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
	plot(0, 0, xlim=c(0, nrow(out)), ylim=c(0,1), xaxs="i", yaxs="i", xaxt="n", xlab="", ylab=ylab, type="n", cex.lab=cex.lab, font.lab=font.lab)
	for (i in seq(N)) {
		x <- rep(i, 5) - c(1,0,0,1,1)
		for (j in seq(K)) {
			xy <- cbind(x, csum[i,j+c(0,0,1,1,0)])
			polygon(xy, border=color[i,j], col=color[i,j])
		}
	}
	if (isTRUE(show.ind.label)) {
		mtext(text=rownames(out)[ord], side=1, line=0.25, at=seq(0.5, nrow(out)-0.5, by=1), cex=cex.axis, las=2)
	}
	if (!missing(expr)) {
		if (is.expression(expr)) {
			eval(expr)
		}
	}
	if (isTRUE(device == "pdf")) {
		dev.off()
	}
}



### opacity ###
# based on add_trans at: https://rdrr.io/github/mcsiple/mmrefpoints/man/add_trans.html
# Arguments:
# color - input color
# alpha - degree of opacity (from 0 to 1)	
# Details:
# this function defines opacity of a color, 0 being fully transparent and 1 being fully visable
# works with either color and opacity a vector of equal length, or one of the two of length 1.

opacity <- function(color, alpha) {
	num2hex <- function(x) {
		hex <- unlist(strsplit("0123456789ABCDEF",split=""))
		return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
	}
	if (length(color) != length(alpha) & !any(c(length(color), length(alpha)) == 1)) stop("Vector lengths not correct")
	if (length(color) == 1 & length(alpha) > 1) color <- rep(color, length(alpha))
	if (length(alpha) == 1 & length(color) > 1) alpha <- rep(alpha, length(color))
	alpha <- as.integer(255 * alpha)
	rgb <- rbind(grDevices::col2rgb(color), alpha)
	res <- paste("#", apply(apply(rgb, 2, num2hex), 2, paste, collapse=""), sep="")
	return(res)	
}


custom_str_palette <- function(str_res, info, palette, ind="ID", pop="Species") {
	trubetskoy <- c(Red="#e6194b", Green="#3cb44b", Yellow="#ffe119", Blue="#0082c8", Orange="#f58231", Purple="#911eb4", Cyan="#46f0f0", Magenta="#f032e6", Lime="#d2f53c", Pink="#fabebe", Teal="#008080", Lavender="#e6beff", Brown="#aa6e28", Beige="#fffac8", Maroon="#800000", Mint="#aaffc3", Olive="#808000", Coral="#ffd8b1", Navy="#000080", Grey="#808080", White="#FFFFFF", Black="#000000")
	K <- attr(str_res,"K")
	cluster <- apply(str_res[,seq(K)], 1, which.max)
	population <- info[match(names(cluster), info[,ind]), pop]
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


