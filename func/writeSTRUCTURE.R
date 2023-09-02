read_str_data <- function(x, row.names=1) {
	str <- read.delim(x, skip=1, header=FALSE, row.names=row.names)
	colnams <- paste(rep(unlist(strsplit(readLines(x, n=1), "[[:blank:]]+")), each=2), c(0,1), sep="_")
	add <- ncol(str) - length(colnams)
	if (add > 0) {
		atr <- str[,seq(add)]
		str <- str[,-seq(add)]
	} else {
		atr <- NULL
	}
	colnames(str) <- colnams
	#return(str[,grepl("^[-]*[[:digit:]]", str[1,])])
	attr(str, "info") <- atr
	return(str)
}


# it supposes phased diploid genotypes labelled by individual name + allele identifier specified by the argument 'allele'

write_str_data <- function(loci, file=NULL, return=FALSE, allele="_[[:digit:]]+$", linkage=FALSE, PopData=NULL, PopFlag=1, structure=NULL, toupper=FALSE) {
if (is.null(structure)) {
	if (inherits(loci, "matrix")) {
		locnam <- colnames(loci)
		if (is.null(locnam)) {
			locnam <- paste0("L", formatC(seq(ncol(loci)), format="d", flag=0, width=nchar(ncol(loci))))
		}
		loci <- setNames(lapply(seq(ncol(loci)), function(i) loci[,i,drop=FALSE]), locnam)
	} else {
		if (is.null(names(loci))) {
			names(loci) <- paste0("L", formatC(seq_along(loci), format="d", flag=0, width=nchar(length(loci))))
		}
		for (i in seq_along(loci)) {
			if (is.null(colnames(loci[[i]])) & ncol(loci[[i]]) > 1) {
				colnames(loci[[i]]) <- paste0(names(loci)[i], "S", formatC(seq(ncol(loci[[i]])), format="d", flag=0, width=nchar(ncol(loci[[i]]))))
			} else if (is.null(colnames(loci[[i]])) & ncol(loci[[i]]) == 1) {
				colnames(loci[[i]]) <- names(loci)[i]
			}
		}
	}
	char <- Filter(function(x) sign(x) == 1, suppressWarnings(as.numeric(loci[[1]])))
	type <- ifelse(length(char) == 0, "snps", "num")
	if (isTRUE(toupper) | inherits(loci[[1]], "DNAbin")) {
		loci <- lapply(loci, toupper)
	}
	for (i in seq_along(loci)) {
		loci[[i]] <- loci[[i]][order(rownames(loci[[i]])),,drop=FALSE]
		if (type == "snps") {
			for (j in seq(ncol(loci[[i]]))) {
				loci[[i]][,j] <- match(loci[[i]][,j], c("A","C","G","T"))
			}
		}
		ind <- gsub(allele, "", rownames(loci[[i]]))
		loc <- paste(rep(names(loci[i]), each=2), c(0,1), sep="_")
		loci[[i]] <- matrix(unlist(split(loci[[i]], ind)), length(ind)/2, length(loc), byrow=TRUE, dimnames=list(unique(ind), loc))
	}
	rows <- sort(unique(unlist(lapply(loci, rownames))))
	structure <- vector("list", length(loci))
	for (i in seq_along(structure)) {
		structure[[i]] <- matrix(-9, length(rows), ncol(loci[[i]]))
		structure[[i]] <- loci[[i]][match(rows, rownames(loci[[i]])),]
		structure[[i]][,] <- match(structure[[i]], sort(na.omit(unique(as.vector(structure[[i]])))))
		dimnames(structure[[i]]) <- list(rows, colnames(loci[[i]]))
	}
	if (isTRUE(linkage)) {
		linkage <- paste(unlist(lapply(structure, function(x) rep(c(-1,0), c(1, ncol(x)/2 - 1)))),  collapse="\t")
	}
	structure <- do.call(cbind, structure)
	structure[is.na(structure)] <- "-9"
}
	if (isFALSE(linkage)) {
		linkage <- NULL
	} else {
		linkage <- paste(linkage,  collapse="\t")
	}
	if (!is.null(PopData)) {
		matching <- match(rownames(structure), PopData[,1])
		PopData <- PopData[matching,2]
		PopFlag <- rep_len(as.numeric(PopFlag), length(PopData))[matching]
	} else {
		PopFlag <- NULL
	}
	locnames <- paste(unique(gsub("_[[:digit:]]+$", "", colnames(structure))),  collapse="\t")
	header <- cbind(rownames(structure), PopData, PopFlag)
	input <- cbind(header, structure)
	input <- apply(input, 1, paste,  collapse="\t")
	input <- c(locnames, linkage, input)
	if (!is.null(file)) {
		writeLines(input, con=file)
	}
	if (isTRUE(return)) {
		return(structure)
	}
}






parse_str_params <- function(df, templates =c("mainparams", "extraparams"), outputs=c("main.txt", "extra.txt")) {
	names(df) <- c("Arg", "Value")
	df$Arg <- as.character(paste0("#define[[:blank:]]+", toupper(df$Arg)))

	main <- readLines(templates[1])
	hits <- lapply(df$Arg, grep, x=main)
	mainpar <- which(sapply(hits, length) == 1)
	for (i in mainpar) {
		legend <- sub("^#[[:alnum:][:blank:]]+\\/", "/", main[hits[[i]]])
		main[hits[[i]]] <- paste(sub("[[:blank:]]+", " ", df$Arg[i], fixed=TRUE), "   ", df$Value[i], "   ", legend)
	}
	writeLines(main, con=outputs[1])
	df <- df[-mainpar,,drop=FALSE]
	extra <- readLines(templates[2])
	if (nrow(df) > 0) {
		hits <- lapply(df$Arg, grep, x=extra)
		extrapar <- which(sapply(hits, length) == 1)
		for (i in extrapar) {
			legend <- sub("^#[[:alnum:][:blank:]]+\\/", "/", extra[hits[[i]]])
			extra[hits[[i]]] <- paste(sub("[[:blank:]]+", " ", df$Arg[i], fixed=TRUE), "   ", df$Value[i], "   ", legend)
		}
	}
	writeLines(extra, con=outputs[2])	
}


### read_str_output ###
# file: name of the file with output of STRUCTURE (so far it works for admixture model with no POPINFO prior)

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


writeGeno <- function(loci, file=NULL, return=FALSE) {
    alleles <- vector("list", length(loci))
    for (i in seq_along(alleles)) {
        bp <- matrix(toupper(loci[[i]]) %in% c("A","C","G","T"), nrow(loci[[i]]), ncol(loci[[i]]))
        sq <- loci[[i]][,apply(bp, 2, all), drop=FALSE]
        sq <- apply(sq, 1, paste, collapse="")
        sq <- sq[order(rownames(loci[[i]]))]
        ind <- gsub("_[[:digit:]]+$", "", rownames(loci[[i]]))
        alleles[[i]] <- setNames(match(sq, unique(sq)), ind)
    }
    rows <- sort(unique(unlist(lapply(alleles, names))))
    geno <- vector("list", length(alleles))
    for (i in seq_along(alleles)) {
        tab <- table(alleles[[i]], names(alleles[[i]]))
        geno[[i]] <- matrix(9, nrow(tab), length(rows), dimnames=list(NULL, rows))
        geno[[i]][, match(colnames(tab), rows)] <- tab
    }
    geno <- do.call(rbind, geno)
    geno <- apply(geno, 1, paste, collapse="")
    if (!is.null(file)) {
        writeLines(geno, con=file)
    }
    if (isTRUE(return)) {
        return(geno)
    }
}

