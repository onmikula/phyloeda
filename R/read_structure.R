#' Reading STRUCTURE data file.
#' 
#' @description
#' Reads STRUCTURE input data file.
#'
#' @param file character string, name of the data file.
#' @param label logical, whether the first column contains individual labels, default is `TRUE`.
#' @param popdata,popflag,locdata logical, whether the initial columns contain information on population of origin (`popdata`),
#'   whether to use it in clustering (`popflag`) and information on "locality", i.e. the clustering affecting prior (`locdata`).
#' @param markernames logical, whether the first row contains marker (=locus) names, default is `TRUE`.
#' @param recessive,linkage logical, whether the initial rows contain information which allele is recessive (`recessive`),
#'   and what are genomic distances separating the markers (`linkage `).
#' @param phase logical, whether rows with individual genotypes are alternated by rows with phase information.
#' @param ploidy numeric, the ploidy, i.e., the number of rows per individual genotype, default is 2.
#' @param allele optional character vector with identifiers used to distinguish individual alleles, e.g., `c("_0", "_1")`,
#'   used only if there is any phase information.
#' @returns A matrix with genotypes, associated with an attributes `pops` (population information from the first columns),
#'   `markers` (information about markers the first rows) and `phase` with information about phasing of individual genotypes.
#' @details The way STRUCTURE data files are formatted requires explicit specification, which rows and columns are included.
#'   By default, the function assumes the file to contain individual and marker names (`label=TRUE` and `markernames=TRUE`),
#'   but no other information. The function can handle data of arbitrary ploidy, if it does not vary across individuals.
#'   Currently, it does not support formatting with one row per individual diploid genotype, although this can be by-passed
#'   by setting `ploidy=1`.
#' @export

read_str_data <- function(file, label=TRUE, popdata=FALSE, popflag=FALSE, locdata=FALSE, markernames=TRUE, recessive=FALSE, linkage=FALSE, phase=FALSE, ploidy=2, allele) {
	str <- readLines(file)
	str <- gsub("^[[:blank:]]+|[[:blank:]]+$", "", str)
	str <- str[nchar(str) > 0]
	str <- strsplit(str, "[[:blank:]]+")
	if (isTRUE(markernames)) { colnams <- str[1]; str <- str[-1] }
	if (isTRUE(recessive)) { rec <- as.integer(str[1]); str <- str[-1] } else { rec <- NULL }
	if (isTRUE(linkage)) { dst <- as.integer(str[1]); str <- str[-1] } else { dst <- NULL }
	if (isTRUE(phase)) {
		phaserows <- seq_along(str) %% (ploidy + 1) == 0
		phase <- do.call(rbind, str[phaserows])
		str <- str[-phaserows]
	} else {
		phase <- NULL
	}
	str <- do.call(rbind, str)
	if (isTRUE(markernames)) { colnames(str) <- colnams }
	if (isTRUE(label)) { rownames(str) <- str[,1]; str <- str[,-1] }
	ind <- seq(nrow(str)) %% ploidy == 0
	if (isTRUE(popdata)) { pop <- str[ind,1]; str <- str[,-1] } else { pop <- NULL }
	if (isTRUE(popflag)) { usepop <- str[ind,1]; str <- str[,-1] } else { usepop <- NULL }
	if (isTRUE(locdata)) { loc <- str[ind,1]; str <- str[,-1] } else { loc <- NULL }
	pops <- cbind(popdata=pop, popflag=usepop, locdata=loc)
	markers <- cbind(recessive=rec, linkage=dst)
	if (!missing(allele) & !is.null(phase)) {
		rownames(str) <- paste0(rownames(str), allele[seq(ploidy)])
	}
	mode(str) <- "integer"
	attr(str, "pops") <- pops
	attr(str, "markers") <- markers
	attr(str, "phase") <- phase
	return(str)
}



#' Reading STRUCTURE output file.
#' 
#' @description
#' Reads STRUCTURE output data file.
#'
#' @param file character string, name of the output file.
#' @returns A matrix with estimated ancestry proportions and their credibility intervals.
#' @details So far it works only for admixture model with no POPINFO prior.
#' @export

read_str_output <- function(file) {
	x <- suppressWarnings(readLines(file))
	x <- x[nchar(x) > 0]
	K <- as.numeric(sub("MAXPOPS=", "", regmatches(x, regexpr("MAXPOPS=[[:digit:]]+", x))))
	nind <- as.numeric(sub("NUMINDS=", "", regmatches(x, regexpr("NUMINDS=[[:digit:]]+", x))))
	start <- grep("Inferred ancestry of individuals", x) + 2
	end <- start + nind - 1
	qmat <- x[start:end]
	qmat <- gsub("^[[:blank:]]+", "", qmat)
	qmat <- do.call(rbind, strsplit(qmat, "[[:blank:]]+"))
	div <- which(qmat[1,] == ":")
	qdf <- data.frame(Label=qmat[,2], Miss=as.numeric(gsub("[\\(\\)]", "", qmat[,3])))
	est <- qmat[,div+(1:K)]
	mode(est) <- "numeric"
	rownames(est) <- qdf$Label
	colnames(est) <- paste0("Cluster", formatC(seq(K), format="d", flag=0, width=nchar(K)))
	for (i in seq(K)) {
		int <- gsub("[\\(\\)]", "", qmat[,K + div + i])
		int <- do.call(rbind, lapply(strsplit(int, ","), as.numeric))
		colnames(int) <- paste(colnames(est)[i], c("low","upp"), sep="_")
		est <- cbind(est, int)
	}
	est <- as.data.frame(est)
	attr(est, "K") <- K
	attr(est, "Missing") <- qdf$Miss
	return(est)
}

