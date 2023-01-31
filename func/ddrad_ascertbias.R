# Function:
#	ddrad_ascertbias: calculates statistics necessary for the analysis of ascertainment bias across a ddRADseq loci
# Arguments:
#	loci: list of loci or snps
#	pop: optional, data frame with two columns: the first containing individual IDs without allelic specification, the second classification to a taxon (population, species etc.)
#	thinning - consider only every x-th locus (thinning = x)
# Value:
#	a list with components occupancy (% of successfully sequenced individuals), mdist (mean uncorrected genetic distance), nalleles (the number of distinct haplotypes), nsnps (the number of SNPs), ntaxa (the number of taxa represented by the locus), mintra (mean intra-population / intra-specific uncorrected genetic distance), thinning (the thinning applied). The components with statistics are vectors of the same length as the number of loci after thinning.

ddrad_ascertbias <- function(loci, pop, thinning=1) {
	loci <- loci[seq_along(loci) %% thinning == 0]
    if (class(loci[[1]])[1] == "matrix") {
        loci <- lapply(loci, ape::as.DNAbin)
    }
	seqnam <- sort(unique(unname(unlist(lapply(loci, rownames)))))
	popanalysis <- !missing(pop)
	occupancy <- sapply(loci, nrow) / length(seqnam)
	nalleles <- setNames(rep(1, length(loci)), names(loci))
	nsnps <- matrix(0, length(loci), 2, dimnames=list(names(loci), c("snp", "pis")))
	mdist <- setNames(rep(NA, length(loci)), names(loci))
    if (isTRUE(popanalysis)) {
        ntaxa <- setNames(rep(NA, length(loci)), names(loci))
        mintra <- setNames(rep(NA, length(loci)), names(loci))
        indnam <- setNames(lapply(pop[,1], grep, x=seqnam), pop[,1])
        indnam <- rep(names(indnam), sapply(indnam, length))[order(unlist(indnam))]
    } else {
        ntaxa <- NULL
        mintra <- NULL
    }
	for (i in seq_along(loci)) {
		dst <- ape::dist.dna(loci[[i]], model="raw", pairwise.deletion=TRUE)
		mdist[i] <- mean(dst)
		snps <- find_snps(loci[[i]])
		nsnps[i,"snp"] <- length(snps$snp)
		if (nsnps[i,"snp"] > 0) {
			nsnps[i,"pis"] <- sum(snps$pis)
			un <- unique(toupper(loci[[i]][,snps$snp,drop=FALSE]))
			un[!un %in% c("A","C","G","T")] <- "N"
			un <- un[order(rowSums(un == "N")),,drop=FALSE]		
			rows <- unname(apply(un != "N", 1, all))
			for (j in which(!rows)) {
				subs <- un[c(j, which(rows)),,drop=FALSE]
				subs <- subs[,colSums(subs == "N") == 0,drop=FALSE]
				rows[j] <- !duplicated(subs, fromLast=TRUE)[1]
			}
			nalleles[i] <- sum(rows, na.rm=TRUE)
		}
        if (isTRUE(popanalysis)) {
            ind <- indnam[match(rownames(loci[[i]]), seqnam)]
            popnam <- pop[match(ind, pop[,1]),2]
            ntaxa[i] <- length(unique(popnam))
            intra <- !diag(length(ind)) & outer(popnam, popnam, "==")
            if (any(intra == TRUE)) {
                mintra[i] <- mean(dst[intra], na.rm=TRUE)
            }
        }
	}
	result <- list(occupancy=occupancy, mdist=mdist, nalleles=nalleles, nsnps=nsnps, ntaxa=ntaxa, mintra=mintra, thinning=thinning)
	return(result)
}



# Function:
#	plotAscertBias: calculates statistics necessary for the analysis of ascertainment bias across a ddRADseq loci
# Arguments:
#	x: name of x-axis variable, e.g. 'occupancy'
#	y: name of y-axis variable, e.g. 'mdist'
#	ascert: output of ddrad_ascertbias
#	cumulative: if TRUE (default), the plot shows statistics for increasingly filtered sets of loci
#	xlab, ylab, main: labels of axes and the whole plot
#	cex, cex.lab, cex.axis: size of points, axis labels and tick labels
#	return: whether to return a list with binned statistics ('classes'), their medians ('medians'), break points between bins ('breaks') and the numbers of loci contained in the bins ('nloci')

plotAscertBias <- function(x, y, ascert, cumulative=TRUE, xlab="", ylab="", main="", cex=1.5, cex.lab=1.5, cex.axis=1.25, return=FALSE) {
	makebreaks <- function(ascert, what) {
		breaks <- c(seq(0, 1, by=0.05), 1.0000001)
		bins <- .bincode(ascert[[what]], breaks, right=FALSE)
		binlab <- sort(unique(bins))
		breaks <- breaks[binlab]
		names(breaks) <- formatC(breaks, format="f", digit=2)
		return(list(breaks=breaks, bins=bins, n=length(breaks)))
	}
	bpts <- makebreaks(ascert, what=x)
	classes <- split(ascert[[y]], bpts$bins)
	if (isTRUE(cumulative)) {
		for (i in (length(classes)-1):1) {
			classes[[i]] <- c(classes[[i]], classes[[i+1]])
		}
	}
	medians <- sapply(classes, median)
	nloci <- round(ascert$thinning * sapply(classes, length))
	xlim <- c(0.5, bpts$n+0.5)
	ylim <- range(ascert[[y]], na.rm=TRUE)
	par(mar=c(5.1,4.1,3.1,4.1))
	plot(seq(bpts$n), medians, xlim=xlim, ylim=ylim, xlab="", ylab="", main=main, type="n", xaxt="n", cex.axis=cex.axis)
	axis(1, at=seq(bpts$n), labels=names(bpts$breaks), cex.axis=cex.axis)
	mtext(xlab, side=1, line=3, cex=cex.lab)
	mtext(ylab, side=2, line=2.25, cex=cex.lab)
	mtext("no. of loci", side=4, line=2.25, cex=cex.lab)
	for (i in seq(bpts$n)) {
		tryCatch(vioplot::vioplot(classes[[i]], add=TRUE, pchMed=NA, at=i), error=function(e) NULL)
	}
	points(seq(bpts$n), medians, pch=21, bg="white", col="black", cex=cex)
	zlim <- range(nloci)
	ratio <- diff(ylim) / diff(zlim)
	rescaled <- ratio * (nloci - zlim[1]) + ylim[1]
	at <- axTicks(4, c(zlim, 5))
	axis(4, at=ratio * (at - zlim[1]) + ylim[1], labels=at, las=3, cex.axis=0.9*cex.axis)
	points(seq(bpts$n), rescaled, pch=21, col=3, bg=3, cex=cex)
	if (isTRUE(return)) {
		return(list(medians=medians, classes=classes, breaks=bpts, nloci=nloci))
	}
}
