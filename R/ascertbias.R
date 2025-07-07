#' Ascertainment bias.
#' 
#' @description
#' Calculates statistics helping to detect of ascertainment bias due to different occupancy across loci.
#'
#' @param loci a list of locus-specific sequence alignments of class `matrix` or `DNAbin`.
#' @param info a data frame or matrix listing individual labels (1st column) and group labels (2nd column).
#' @param allele regular expression, allele identifier in rownames of locus-specific alignments.
#' @param thin a proportion of loci to be retained or a subsampling frequency (if `thin > 1`).
#' @returns A list with components `occupancy` (proportion of successfully sequenced individuals),
#'   `mdist` (mean uncorrected genetic distance), `nalleles` (the number of distinct haplotypes),
#'   `nsnps` (the number of SNPs), `npis` (the number of parsimony-informative SNPs),
#'   `heterozygosity` (mean heterozygosity), `spoccupancy` (proportion of species or other taxa represented
#'   at the locus), `mintra` (mean intra-specific uncorrected genetic distance) and `thin` (the `thin` argument).
#'   The components are numeric vectors of the same length as the number of loci after thinning.
#' @export

ascertbias <- function(loci, info, allele="_[012]$", thin=1) {

	loci <- subset_loci(loci, sapply(loci, nrow) > 0)
	if (all(nchar(names(loci)) == 0)) {
		names(loci) <- paste0("L", formatC(seq_along(loci), format="d", flag=0, width=nchar(length(loci))))
	}
	if (thin != 1) {
		loci <- thinning(loci, freq=thin)
	}
	snps <- find_snps(loci, allele=allele)
    if (inherits(loci[[1]], "matrix")) {
        loci <- lapply(loci, ape::as.DNAbin)
    }

	occ <- occupancy(loci)
	dst <- lapply(loci, ape::dist.dna, model="raw", pairwise.deletion=TRUE)
	nalleles <- setNames(rep(1, length(loci)), names(loci))
	nsnps <- matrix(0, length(loci), 2, dimnames=list(names(loci), c("snp", "pis")))
	mdist <- setNames(rep(NA, length(loci)), names(loci))
	hetzyg <- setNames(rep(0, length(loci)), names(loci))
	for (i in seq_along(loci)) {
		loci[[i]] <- loci[[i]][order(rownames(loci[[i]])),,drop=FALSE]
		mdist[i] <- mean(dst[[i]], na.rm=TRUE)
		nsnps[i,"snp"] <- length(snps[[i]]$snp)
		if (nsnps[i,"snp"] > 0) {
			varsites <- toupper(loci[[i]])[,snps[[i]]$snp,drop=FALSE]
			varsites[!varsites %in% c("A","C","G","T")] <- "N"
			nsnps[i,"pis"] <- sum(snps[[i]]$pis)
			un <- unique(varsites)
			un <- un[order(rowSums(un == "N")),,drop=FALSE]		
			rows <- unname(apply(un != "N", 1, all))
			for (j in which(!rows)) {
				subs <- un[c(j, which(rows)),,drop=FALSE]
				subs <- subs[,colSums(subs == "N") == 0,drop=FALSE]
				rows[j] <- !duplicated(subs, fromLast=TRUE)[1]
			}
			nalleles[i] <- sum(rows, na.rm=TRUE)
			varsites <- varsites[grepl(allele, rownames(varsites)),,drop=FALSE]
			if (nrow(varsites) >= 2) {
				odd <- seq(nrow(varsites)) %% 2 == 1
				halfmissing <- sum(((varsites[odd,] == "N") + (varsites[!odd,] == "N")) == 1)
				heterozygose <- sum(varsites[odd,] != varsites[!odd,])
				odd <- seq(nrow(loci[[i]])) %% 2 == 1
				total <- sum(loci[[i]][odd,] != "N" & loci[[i]][!odd,] != "N")
				hetzyg[i] <- (heterozygose - halfmissing) / total
				if (hetzyg[i] == Inf) {
					hetzyg[i] <- NA
				}
			} else {
				hetzyg[i] <- NA
			}
		}
 	}

	group <- !missing(info)
	if (isTRUE(group)) {
		grocc <- spoccupancy(loci, info, allele)
		seqnam <- sort(unique(unname(unlist(lapply(loci, rownames)))))
		groups <- info[match(sub(allele, "", seqnam), info[,1]),2]
		dst <- lapply(dst, as.matrix)
        mintra <- setNames(rep(NA, length(loci)), names(loci))
 		for (i in seq_along(loci)) {
 			grp <- groups[match(rownames(dst[[i]]), seqnam)]
            intra <- !diag(length(grp)) & outer(grp, grp, "==")
            if (any(intra == TRUE)) {
                mintra[i] <- mean(dst[[i]][intra], na.rm=TRUE)
            } 			
 		}
 	} else {
		grocc <- NULL
        mintra <- NULL
	}
	
	result <- list(occupancy=occ, mdist=mdist, nalleles=nalleles, nsnps=nsnps[,"snp"], npis=nsnps[,"pis"], heterozygosity=hetzyg, spoccupancy=grocc, mintra=mintra, thin=thin)
	return(result)
	
}



#' Ascertainment bias plot.
#' 
#' @description
#' Plots results of the ascertainment bias analysis.
#'
#' @param var a name of indicator (y-axis) variable; must be a name of component of the list returned by [ascertbias],
#'   i.e. one of the following: 'mdist', 'nalleles', 'nsnps', 'npis', 'heterozygosity', 'spoccupancy', 'mintra'
#'   or (possibly) 'occupancy', or their unambiguous abbreviation.
#' @param ascert output of [ascertbias].
#' @param bin a name of binning (x-axis) variable. Options are the same as in `var`, but the default is 'occupancy'.
#' @param cumulative logical, whether to plot statistics for increasingly filtered sets of loci
#'   rather than for different occupancy bins. The default is `TRUE`.
#' @param xlab,ylab,main labels of axes and the whole plot; if not specified, `xlab` and `ylab` are given
#'   pre-defined values given the specification of `bin` (x-axis) and `var` (y-axis) variables.
#' @param cex,cex.lab,cex.axis size of points, axis labels and tick labels.
#' @param return logical, whether to return a list with binned statistics (`classes`), their medians (`medians`),
#'   break points between bins (`breaks`) and the numbers of loci contained in the bins (`nloci`).
#' @param device either `"quartz"`, `"x11"` (= `"X11"`) or `"pdf"`.
#' @param file character strigm the name of pdf file if `device == "pdf"`.
#' @export

plot_ascertbias <- function(var, ascert, bin="occupancy", cumulative=TRUE, xlab, ylab, main="", cex=1.5, cex.lab=1.5, cex.axis=1.25, return=FALSE, device, file) {
	
	makebreaks <- function(ascert, what) {
		breaks <- c(seq(0, 1, by=0.05), 1.0000001)
		bins <- .bincode(ascert[[what]], breaks, right=FALSE)
		binlab <- sort(unique(bins))
		breaks <- breaks[binlab]
		names(breaks) <- formatC(breaks, format="f", digit=2)
		return(list(breaks=breaks, bins=bins, n=length(breaks)))
	}
	
	var <- base::pmatch(var, names(ascert))
	if (is.na(var)) stop("'var' is not a name of 'ascert' component or its unambiguous abbreviation")
	bin <- base::pmatch(bin, names(ascert))
	if (is.na(bin)) stop("'bin' is not a name of 'ascert' component or its unambiguous abbreviation")
	bpts <- makebreaks(ascert, what=bin)
	classes <- split(ascert[[var]], bpts$bins)
	if (isTRUE(cumulative)) {
		for (i in (length(classes) - 1):1) {
			classes[[i]] <- c(classes[[i]], classes[[i + 1]])
		}
	}
	medians <- sapply(classes, median, na.rm=TRUE)
	nloci <- round(ascert$thin * sapply(classes, length))
	xlim <- c(0.5, bpts$n + 0.5)
	ylim <- range(ascert[[var]], na.rm=TRUE)
	ylabs <- c(occupancy="occupancy", mdist="mean pairwise distance", nalleles="no. of alleles", nsnps="no. of SNPs", npis="no. of parsimony-informative SNPs", heterozygosity="mean heterozygosity", spoccupancy="species-level occupancy", mintra="mean intra-specific pairwise distance")
	if (missing(ylab)) ylab <- ylabs[var]
	if (missing(xlab)) xlab <- paste("min.", ylabs)[bin]

	# plotting
	if (missing(device)) {
		device <- ifelse(.Platform$OS.type == "unix", "quartz", "x11")
	} else if (!is.null(device)) {
		if (isTRUE(device == "pdf")) {
			if (missing(file)) {
				file <- "ascertbias.pdf"
			}
			pdf(file, width=7, height=7)
		} else {
			match.fun(tolower(device))(width=7, height=7)	
		}		
	}

	par(mar=c(5.1,4.1,3.1,4.1))
	plot(seq(bpts$n), medians, xlim=xlim, ylim=ylim, xlab="", ylab="", main=main, type="n", xaxt="n", cex.axis=cex.axis)
	axis(1, at=seq(bpts$n), labels=names(bpts$breaks), cex.axis=cex.axis, las=2)
	mtext(xlab, side=1, line=3.6, cex=cex.lab)
	mtext(ylab, side=2, line=2.35, cex=cex.lab)
	mtext("no. of loci", side=4, line=2.35, cex=cex.lab)
	for (i in seq(bpts$n)) {
		tryCatch(vioplot::vioplot(classes[[i]], add=TRUE, pchMed=NA, at=i), error=function(e) NULL)
	}
	points(seq(bpts$n), medians, pch=21, bg="white", col="black", cex=cex)
	zlim <- range(nloci)
	ratio <- diff(ylim) / diff(zlim)
	rescaled <- ratio * (nloci - zlim[1]) + ylim[1]
	at <- as.integer(axTicks(4, c(zlim, 5)))
	axis(4, at=ratio * (at - zlim[1]) + ylim[1], labels=at, las=3, cex.axis=0.9*cex.axis)
	points(seq(bpts$n), rescaled, pch=21, col=3, bg=3, cex=cex)
	if (isTRUE(return)) {
		return(list(medians=medians, classes=classes, breaks=bpts, nloci=nloci))
	}

	if (isTRUE(device == "pdf")) {
		invisible(dev.off())
	}

}
