#' Writing STRUCTURE data file.
#' 
#' @description
#' Writes STRUCTURE data file.
#'
#' @param loci input data, a matrix or a list of multiple sequence alignments.
#' @param file character string, name of the output file.
#' @param popdata a data frame or a matrix with classification of individuals (1st column) into populations (2nd column).
#' @param popflag logical, whether to use the popdata, can be a vector of the same length as the number of individuals.
#' @param locdata a data frame or a matrix with classification of individuals (1st column) into a priori clusters (2nd column).
#' @param recessive logical, whether to specify recessive alleles, or integer, giving the recessive alleles for each locus.
#' @param linkage logical, whether to apply linkage model, or numeric, giving distances between successive loci.
#'   If distances are specified, the first locus in each linkage group has value -1.
#' @param haploid a list of two vectors containing individual names (1st vector) and loci (2nd vector) with haploid genotypes.
#' @param allele regular expression with allele identifier.
#' @param return logical, whether to return genotype matrix
#' @returns A matrix with genotypes and linkage information coded as in STRUCTURE input file.
#' @export

write_str_data <- function(loci, file=NULL, popdata=NULL, popflag=1, locdata=NULL, recessive=FALSE, linkage=FALSE, haploid=NULL, allele="_[012]$", return=FALSE) {

	if (inherits(loci, "matrix")) {
		locnam <- colnames(loci)
		if (is.null(locnam)) {
			locnam <- paste0("L", formatC(seq(ncol(loci)), format="d", flag=0, width=nchar(ncol(loci))))
		}
		# tady pridat moznost partitioning podle argumentu 'part'
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
		if (isFALSE(linkage) & any(sapply(loci,ncol) > 1)) {
			linkage <- TRUE
		}
	}

	nloci <- length(loci)
	char <- Filter(function(x) sign(x) == 1, suppressWarnings(as.numeric(loci[[1]])))
	type <- ifelse(length(char) == 0, "snps", "num")
	if (type == "snps") {
		for (i in 1:nloci) {
			loci[[i]] <- toupper(loci[[i]])
			for (j in seq(ncol(loci[[i]]))) {
				loci[[i]][,j] <- match(loci[[i]][,j], c("A","C","G","T"))
			}
		}
	}
	ind <- lapply(lapply(loci, rownames), function(x) sub(allele, "", x))
	for (i in 1:nloci) {
		loci[[i]] <- loci[[i]][order(ind[[i]]),,drop=FALSE]
		ind[[i]] <- sort(ind[[i]])
	}
	indnam <- sort(unique(unlist(ind)))
	indploidy <- lapply(ind, table)
	maxploidy <- setNames(numeric(length(indnam)), indnam)
	for (i in seq_along(indploidy)) {
		maxploidy[names(indploidy[[i]])] <- mapply(max, maxploidy[names(indploidy[[i]])], indploidy[[i]])
	}
	part <- cumsum(sapply(loci, ncol))
	part <- cbind(c(1, part[-length(part)]+1), part)
	rownam <- paste(rep(names(maxploidy), maxploidy), unlist(lapply(maxploidy, seq)), sep="_")
	structure <- matrix(NA, length(rownam), max(part), dimnames=list(rownam, unlist(lapply(loci, colnames)))) 
	for (i in 1:nloci) {
		irownam <- paste(rep(names(indploidy[[i]]), indploidy[[i]]), unlist(lapply(indploidy[[i]], seq)), sep="_")
		structure[rownam %in% irownam, part[i,1]:part[i,2]] <- loci[[i]]
	}

	if (!is.null(haploid)) {
		haploid <- haploid[order(sapply(haploid, "[[", 1) %in% indnam, decreasing=TRUE)]
		haploid <- setNames(haploid, c("ind", "loc"))
		haploid$ind <- intersect(indnam, haploid$ind)
		haploid$loc <- intersect(colnames(structure), haploid$ind)
		locnam <- paste0("^", gsub("\\|", "\\|", gsub("\\.", "\\.", haploid$loc)))
		j <- sort(unlist(lapply(locnam, grep, x=colnames(structure))))
		if (any(duplicated(j))) {
			stop("locus identifiers in 'haploid' are not unique")
		}
		for (i in seq_along(haploid$ind)) {
			structure[grepl(paste0(haploid$ind[i], "_[2345678]$"), rownam),j] <- NA
		}
		phase <- matrix(0.5, length(indnam), ncol(structure), dimnames=list(indnam, colnames(structure)))
		phase[haploid$ind,j] <- 1
		rownames(phase) <- paste0(indnam, "_PHASE")
		structure <- rbind(structure, phase)
		structure <- structure[order(rownames(structure)),]
		rownames(structure) <- sub("^.+_PHASE", "", rownames(structure))
	}
	structure[is.na(structure)] <- "-9"
	rownames(structure) <- sub("_[[:alnum:]]$", "", rownames(structure))

	if (!is.null(popdata)) {
		matching <- match(rownames(structure), popdata[,1])
		popdata <- popdata[matching,2]
		popflag <- rep_len(as.numeric(popflag), length(popdata))[matching]
		popdata <- replace(popdata, is.na(popdata), "")
		popflag <- replace(popflag, is.na(popflag), "")
	} else {
		popflag <- NULL
	}
	if (!is.null(locdata)) {
		matching <- match(rownames(structure), locdata[,1])
		locdata <- locdata[matching,2]
		locdata <- replace(locdata, is.na(locdata), "")
	}

	locnames <- colnames(structure)
	if (!is.null(locnames)) {
		locnames <- paste(unique(sub("_[[:digit:]]$", "", locnames)),  collapse="\t")
	}
	if (isFALSE(linkage)) {
		linkinfo <- NULL
#		linkinfo <- rep(-1, ncol(structure))
	} else if (isTRUE(linkage)) {
		linkinfo <- unlist(lapply(seq(nloci), function(i) rep(c(-1,0), c(1,diff(part[i,])))))
	} else if (is.numeric(linkage) & is.vector(linkage)) {
		if (length(linkage) == ncol(structure)) {
			linkinfo <- linkage
		} else if (length(linkage) == nloci) {
			linkinfo <- unlist(lapply(seq(nloci), function(i) rep(c(linkage[i],0), c(1,diff(part[i,])))))
		} else {
			warning("linkage information not used as the length of 'linkage' argument does not match the number of markers")
		}	
	}
	if (!is.null(linkinfo)) {
		linkinfo <- paste(linkinfo,  collapse="\t")
	}
	if (isFALSE(recessive)) {
		recinfo <- NULL
	} else {
		recinfo <- paste(recessive,  collapse="\t")
	}
	header <- cbind(rownames(structure), popdata, popflag, locdata)
	header[is.na(header)] <- ""
	input <- cbind(header, structure)
	input <- apply(input, 1, paste,  collapse="\t")
	input <- c(locnames, recinfo, linkinfo, input)

	if (!is.null(file)) {
		writeLines(input, con=file)
	}
	if (isTRUE(return)) {
		return(structure)
	}

}




#' Writing STRUCTURE parameter files.
#' 
#' @description
#' Prepares STRUCTURE parameter files by accommodating their templates.
#'
#' @param df a data frame with argument names (1st column) and values (2nd column).
#' @param templates character vector, names of the template files.
#' @param outputs character vector, names of the parameter files.
#' @export

write_str_params <- function(df, templates=c("mainparams", "extraparams"), outputs=c("main.txt", "extra.txt")) {
	names(df) <- c("Arg", "Value")
	df$Arg <- as.character(paste0("#define[[:blank:]]+", toupper(df$Arg)))
	linkage <- any(df$Arg == "LINKAGE")
	mapdist <- any(df$Arg == "MAPDISTANCES")
	if (linkage & mapdist) {
		if (df$Value[df$Arg == "LINKAGE"] != df$Value[df$Arg == "MAPDISTANCES"]) {
			stop("values of 'LINKAGE' & 'MAPDISTANCES' must agree")
		}
	} else if (linkage) {
		df <- rbind(df, data.frame(Arg="MAPDISTANCES", Value=df$Value[df$Arg == "LINKAGE"]))
	} else if (mapdist) {
		df <- rbind(df, data.frame(Arg="LINKAGE", Value=df$Value[df$Arg == "MAPDISTANCES"]))
	}

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
