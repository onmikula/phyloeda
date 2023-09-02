### FUNCTION
# even_subsampling
### ARGUMENTS
# bed: tab-delimited file in .bed format
# nloci: no. of loci to be retained
# chromsize: 

even_subsampling <- function(bed, nloci, chromsize) {
	prop <- function(x) x / sum(x)
	chrom <- split(bed, bed[,1])
	if (is.character(chromsize) & length(chromsize) == 1) {
		chrsizes <- read.delim(chromsize, header=FALSE)
		chrsizes <- chrsizes[match(names(chrom), sizes[,1]),]
	} else {
		chrsizes <- chromsize[match(names(chrom), chrsizes[,1]),]
	}
	chrsizes <- sort(setNames(chrsizes[,2], chrsizes[,1]), decreasing=TRUE)
	capacity <- nloci * prop(chrsizes)
	chrom <- chrom[names(capacity)]
	available <- sapply(chrom, nrow)

	dif <- available < capacity
	while (any(dif)) {
		capacity <- ifelse(dif, available, capacity)
		capacity[!dif] <- round((nloci - sum(capacity[dif])) * prop(capacity[!dif]))
		dif <- available < capacity
	}

	if (any(capacity < 1 & available >= 1) & nloci > sum(floor(capacity))) {
		dif <- nloci - sum(floor(capacity))
		add <- which(capacity < 1 & available >= 1)
		capacity[add] <- ifelse(seq_along(add) <= min(c(length(add), dif)), 1, 0)
		if (nloci < sum(capacity)) {
			capacity[capacity > 1] <- floor(capacity[capacity > 1])
		}
	} else if (nloci > sum(floor(capacity))) {
		dif <- nloci - sum(floor(capacity))
		add <- which(available > floor(capacity))
		while (dif > 0 & length(add) > 0) {
			plus <- add[which.max(chrsizes[add] / floor(capacity[add]))]
			capacity[plus] <- floor(capacity[plus] + 1)
			dif <- nloci - sum(floor(capacity))
			add <- which(available > floor(capacity))
		}
		capacity <- floor(capacity)
	}

	chrom <- chrom[names(capacity)]
	for (i in seq_along(capacity)) {
		loci <- chrom[[i]]
		n <- nrow(loci)
		dif <- loci[-1,2] + 1 - loci[-n,3]
		dif <- mapply(min, c(dif[1], dif), c(dif, dif[n-1]))
		while (n > capacity[i]) {
			loci <- loci[-which.min(dif),,drop=FALSE]
			n <- n - 1
			dif <- loci[-1,2] + 1 - loci[-n,3]
			dif <- mapply(min, c(dif[1], dif), c(dif, dif[n-1]))
		}
		chrom[[i]] <- loci
	}

	evensample <- do.call(rbind, chrom)
	row.names(evensample) <- NULL
	return(evensample)
	
}
