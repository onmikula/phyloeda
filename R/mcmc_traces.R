#' Reading MCMC traces.
#' 
#' @description
#' Reads posterior sample of trees or parameters obtained in MCMC-based phylogenetic inference.
#'
#' @param file name of the input file.
#' @param comment.char the commenting character.
#' @param row.names,sep,header optional arguments of [utils::read.table],
#'   relevant if the the imput file is in tabular format.
#' @returns An object of class `multiPhylo` or a data frame, depending on the nature of the analysis.
#' @export

read_trace <- function(file, comment.char="#", row.names=NULL, sep="", header=TRUE) {
	b <- gsub("^[[:blank:]]+|[[:blank:]]+$", "", readLines(file, n=100))
	b <- b[nchar(b) > 0]
	isnexus <- any(grepl("^#*NEXUS", b, ignore.case=TRUE))
	isnewick <- any(grepl("^\\(.+\\);", b))
	if (isnexus | isnewick) {
		x <- read_tree(file)
	} else {
		if (isTRUE(nchar(comment.char) == 1)) {
			x <- read.table(file, comment.char=comment.char, header=header, row.names=row.names, sep=sep)
		} else {
			comment.char <- ifelse(substr(comment.char, 1, 1) != "^", paste0("^", comment.char), comment.char)
			x <- gsub("^[[:blank:]]+|[[:blank:]]+$", "", readLines(file))
			x <- x[!grepl(comment.char, x)]
			x <- read.table(text=x, header=header, row.names=row.names, sep=sep)
		}
	}
	return(x)
}



#' Discarding burn-in.
#' 
#' @description
#' Discards an initial, burn-in, part of the posterior sample.
#'
#' @param x an object of class `multiPhylo` or a data frame.
#' @param relburnin burn-in period specified as a proportion of posterior sample.
#' @param absburnin burn-in period specified as the number of sampled values, overrides `relburnin`.
#' @returns A modified input, with the burn-in part discarded.
#' @export

discard_burnin <- function(x, relburnin, absburnin) {
	istree <- inherits(x, "multiPhylo") | inherits(x[[1]], "phylo")
	isvector <- is.numeric(x) & is.vector(x)
	istable <- is.data.frame(x) | is.matrix(x)
	if (istree | isvector) {
		absburnin <- ifelse(missing(absburnin), ceiling(relburnin * length(x)), absburnin)
		x <- x[-seq(absburnin)]
	} else if (istable) {
		absburnin <- ifelse(missing(absburnin), ceiling(relburnin * nrow(x)), absburnin)
		x <- x[-seq(absburnin),,drop=FALSE]
	} else {
		stop("the input trace must contain either phylogenetic trees in 'phylo' format or sampled values in a vector or a table")
	}
	return(x)
}	


#' Combining posterior samples.
#' 
#' @description
#' Combines posterior samples obtained in independent runs of the same analysis.
#'
#' @param x a list of `multiPhylo` objects or a data frames.
#' @returns An object of class `multiPhylo` or a data frame containing pooled posterior samples.
#' @export

combine_runs <- function(x) {
	istree <- inherits(x[[1]], "multiPhylo") | inherits(x[[1]][[1]], "phylo")
	isvector <- is.numeric(x[[1]]) & is.vector(x[[1]])
	istable <- is.data.frame(x[[1]]) | is.matrix(x[[1]])
	if (istree | isvector) {
		x <- Reduce(c, x)
	} else if (istable) {
		x <- do.call(rbind, x)
	} else {
		stop("the input traces must contain either phylogenetic trees in 'phylo' format or sampled values in a vector or a table")
	}
	return(x)
}


#' Thinning posterior samples.
#' 
#' @description
#' Subsamples posterior samples so every i-th value is retaind.
#'
#' @param x an object of class `multiPhylo` or a data frame.
#' @param freq subsampling frequency or a proportion of samples to be retained.
#' @returns An object of class `multiPhylo` or a data frame containing reduced posterior sample.
#' @export

thinning <- function(x, freq=1) {
	if (freq != 1) {
		thin <- floor(ifelse(freq < 1, 1 / freq, freq))
	}
	istree <- inherits(x, "multiPhylo") | inherits(x[[1]], "phylo")
	isvector <- is.vector(x) & !is.data.frame(x) 
	istable <- is.data.frame(x) | is.matrix(x)
	if (istree | isvector) {
		x <- x[as.integer(seq(0, length(x), by=thin))[-1]]
	} else if (istable) {
		x <- x[as.integer(seq(0, nrow(x), by=thin))[-1],,drop=FALSE]
	} else {
		stop("the input traces must contain either phylogenetic trees in 'phylo' format or sampled values in a vector or a table")
	}
	return(x)	
}

