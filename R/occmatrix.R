#' Two-way occupancy.
#' 
#' @description
#' Calculates simultaneously occupancy of loci and coverage of individuals.
#'
#' @param loci a list of locus-specific sequence alignments of class `matrix` or `DNAbin`.
#' @param allele regular expression, allele identifier in rownames of locus-specific alignments.
#' @returns A matrix with rows corresponding to individuals and columns to subsets of loci increasingly filtered
#'   by their occupancy. It contains individual coverages in these subsets as calculated by [coverage].
#' @export

occmatrix <- function(loci, allele="_[012]$") {
	occ <- occupancy(loci)
	vec <- seq(0, 1, by=0.05)
	int <- findInterval(occ, vec, rightmost.closed=FALSE, left.open=FALSE, all.inside=FALSE)
	tab <- table(int)
	indiv <- sort(unique(unlist(lapply(loci, function(x) unique(sub(allele, "", rownames(x)))))))
	occmat <- do.call(cbind, lapply(split(loci, int), coverage, prop=FALSE, allele=allele, indiv=indiv))
	miss <- setdiff(seq(21), as.numeric(colnames(occmat)))
	if (length(miss) > 0) {
		ord <- as.character(sort(as.numeric(c(miss, colnames(occmat)))))
		occmat <- cbind(matrix(0, nrow(occmat), length(miss), dimnames=list(rownames(occmat), miss)), occmat)[,ord]
		tab <- c(setNames(rep(0, length(miss)), miss), tab)[ord]
	}
	occmat <- t(apply(occmat[,21:1], 1, cumsum))[,21:1]
	occmat <- occmat %*% diag(1 / rev(cumsum(rev(tab))))
	occmat[is.nan(occmat) | is.na(occmat)] <- 0
	colnames(occmat) <- formatC(vec, format="f", digit=2)
	return(occmat)
}



#' Occupancy heatmap.
#' 
#' @description
#' Creates a heatmap showing results of two-way occupancy analysis.
#'
#' @param occmat the output of [occmatrix].
#' @param device either `"quartz"`, `"x11"` (= `"X11"`) or `"pdf"`.
#' @param file character strigm the name of pdf file if `device == "pdf"`.
#' @export

occheatmap <- function(occmat, device, file) {
	palette <- c(rgb(1, 1, 1), rgb(1, seq(1,0,-1/6), 0), rgb(1, 0, seq(1/6,1,1/6)), rgb(seq(1-1/6,0,-1/6), 0, 1), rgb(0.1,0.1,0.5))
	vec <- seq(0, 1, by=0.05)
	N <- nrow(occmat)
	hh <- nrow(occmat) * 4.56 / 21 + 0.82

	if (missing(device)) {
		device <- ifelse(.Platform$OS.type == "unix", "quartz", "x11")
	}
	if (isTRUE(device == "pdf") & missing(file)) {
		file <- "occheatmap.pdf"
	}
	if (isTRUE(device == "pdf")) {
		pdf(file, width=7, height=hh)
	} else {
		match.fun(tolower(device))(width=7, height=hh)	
	}
	layout(matrix(1:2, 1, 2), widths=c(6, 1))

	par(mai=c(0.22, 1.42, 0.62, 0.02))
	graphics::image(seq(21), seq(N), t(occmat)[,N:1], zlim=c(1e-08, 1), col=palette, axes=FALSE, ann=FALSE, xaxs="i", yaxs="i", bty="n")
	mtext(text=colnames(occmat), side=3, line=0.5, at=seq(21), cex=1, las=2)
	mtext(text=rev(rownames(occmat)), side=2, line=0.5, at=seq(N), cex=0.7, las=2)
	box <- cbind(c(0.5, 21.5, 21.5, 0.5), c(0.5, 0.5, N+0.5, N+0.5))
	polygon(box, border="black", lwd=2)

	extramai <-  0.5 * (N - 21) * 4.56 / 21
	mai <- c(0.22, 0.5, 0.62, 0.25) + c(extramai, 0)
	par(mai=mai)
	graphics::image(0:1, seq(21), t(vec), zlim=c(1e-08, 1), col=palette, axes=FALSE, ann=FALSE, xaxs="i", yaxs="i", bty="n")
	mtext(text=colnames(occmat), side=2, line=0.25, at=seq(21), cex=0.7, las=2)
	box <- cbind(c(0, 1, 1, 0), c(0.5, 0.5, 21.5, 21.5))
	polygon(box, border="black", lwd=2)

	if (isTRUE(device == "pdf")) {
		invisible(dev.off())
	}
}

