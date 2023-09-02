### FUNCTION
# plot_gencoord
### ARGUMENTS
# bedtable: tab-delimited file in .bed format
# chromsize: tab-delimited file with scaffold/chromosome sizes (1st col. - scaffold name, 2nd col. - size in bp, without header)
# no: max. no. of the largest scaffolds/chromosomes to be plotted (defaults to 20)
# subset
# order
# col
# lwd
# chr_wd
# chr_col
# chr_bg
# subs_col
# subs_lwd
# cex.lab
# device: "quartz", "x11" (or "X11") or "pdf", if NULL, the objects are plotted into the current device
# file: name of pdf file if device == "pdf"
# width: width of graphical device
# height: height of graphical device
# mai: 
# return

plot_gencoord <- function(bedtable, chromsize, no=20, subset=NULL, order=FALSE, col="blue", lwd=1, chr_wd=0.2, chr_col="black", chr_bg="lightgray", subs_col="red", subs_lwd=2*lwd, cex.lab=1.25, device, width=7, height=NA, file=NULL, mai, return=FALSE) {
	bedtable$mid <- bedtable[,2] + (bedtable[,3] - bedtable[,2]) / 2
	chrom <- split(bedtable, bedtable[,1])
	if (is.character(chromsize) & length(chromsize) == 1) {
		sizes <- read.delim(chromsize, header=FALSE)
	} else {
		sizes <- chromsize
	}
	rows <- match(names(chrom), sizes[,1])
	if (any(is.na(rows))) {
		stop("some scaffolds listed in .bed file are missing in the table of scaffold sizes")
	}
	sizes <- sizes[rows,,drop=FALSE]
	
	if (no < length(chrom)) {
		chromnames <- names(chrom)
		chrom <- chrom[order(sapply(chrom, nrow), decreasing=TRUE)][seq(no)]
		chrom <- rev(chrom[order(match(names(chrom), chromnames))]) 
		sizes <- sizes[match(names(chrom), sizes[,1]),,drop=FALSE]
	} else {
		no <- length(chrom)
	}
	
	if (isTRUE(order)) {
		sizeord <- order(sizes[,2])
		chrom <- chrom[sizeord]
		sizes <- sizes[sizeord,,drop=FALSE]
	}
	
	if (!is.null(subset)) {
		subs <- setNames(vector("list", no), names(chrom))
      	for (i in seq(no)) {
			subs[[i]] <- chrom[[i]][,4] %in% subset
      	}
	}

	r <- no / 20
	if (is.na(height)) {
		height <- width * r
	}

	# plotting
	if (missing(device)) {
		device <- ifelse(.Platform$OS.type == "unix", "quartz", "x11")
	}
	if (isTRUE(device == "pdf") & missing(file)) {
		file <- "genomic_positions.pdf"
	}
	if (isTRUE(device == "pdf")) {
		pdf(file, width=width, height=height)
	} else if (!is.null(device)) {
		match.fun(tolower(device))(width=width, height=height)	
	}
	
#	par(mai=c(1.02, 1.22, 0.22, 0.22))
	if (missing(mai)) {
		mai <- c(0.92, 0.92, 0.12, 0.12)
	}
	par(mai=mai)
	plot(1, 1, type="n", bty="n", xlim=c(0, max(sizes[,2])), ylim=c(1 - 2 * chr_wd, no + 2 * chr_wd), yaxt="n", xlab="position [bp]", ylab="", cex.lab=cex.lab)
	mtext(text=rev(names(chrom)), side=2, line=0.1, outer=FALSE, at=seq(no), adj=NA, padj=NA, cex=0.7*cex.lab, col=1, font=1, las=2)
      	for (i in seq(no)) {
		x <- c(0, sizes[i,2])[c(1,2,2,1)]
		y <- rep(i, 4) + c(-chr_wd, -chr_wd, chr_wd, chr_wd)
		polygon(x, y, border=chr_col, col=chr_bg)
		for (j in seq(nrow(chrom[[i]]))) {
			lines(cbind(chrom[[i]]$mid[j], i + c(-chr_wd, chr_wd)), col=col, lwd=lwd)
		}
		if (!is.null(subset)) {
			for (j in which(subs[[i]])) {
				lines(cbind(chrom[[i]]$mid[j], i + c(-chr_wd, chr_wd)), col=subs_col, lwd=subs_lwd)
			}
		}
	}
	
	if (isTRUE(device == "pdf")) {
		dev.off()
	}

	if (isTRUE(return)) {
		res <- do.call(rbind, chrom)
		rownames(res) <- NULL
		return(res)
	}

}

