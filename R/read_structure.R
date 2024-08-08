#' Reading STRUCTURE data file.
#' 
#' @description
#' Reads STRUCTURE input data file.
#'
#' @param file character string, name of the data file.
#' @param row.names a vector of row names or an index of the column containing row names.
#' @returns A matrix with genotypes, associated with an attribute `info` which contains the first columns with metadata.
#' @details It supposes haploid or diploid genotypes.
#' @export

read_str_data <- function(file, row.names=1) {
	str <- readLines(file)
	colnams <- unlist(strsplit(str[1], "[[:blank:]]+"))
	rownams <- sub("\\s.+$", "", str[-1])
	if(!any(duplicated(rownams[nchar(rownams) > 0]))) {
		colnams <- paste(rep(colnams, each=2), c(0,1), sep="_")
	}
	str <- sub("^.+\\s", "", str[-1])
	str <- as.matrix(read.delim(file, skip=1, header=FALSE, row.names=row.names))
	colnams <- unlist(strsplit(readLines(file, n=1), "[[:blank:]]+"))
	add <- ncol(str) - length(colnams)
	if (add > 0) {
		atr <- str[,seq(add)]
		str <- str[,-seq(add)]
	} else {
		atr <- NULL
	}
	colnames(str) <- colnams
	attr(str, "info") <- atr
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

