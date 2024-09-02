#' Read partition file.
#' 
#' @description
#' Reads nexus or RAxML-formatted partition files, including files created by ModelFinder tool of IQ-TREE software.
#'
#' @param file a name of the partition file.
#' @param param logical, whether to read also estimates of substitution model parameters.
#' @param name character string with locus name, relevant when reading in standard output of IQ-TREE
#'   (rather than a partition file it creates).
#' @returns A data frame with columns for the partition name, start, end, codon position
#'   and possibly for components of nucleotide substitution models and estimates of their parameters.
#'   Every row in the data frame corresponds to a single elementary partition, which can be merged with another one,
#'   as indicated by the partition name.
#' @export

read_part <- function(file, param=FALSE, name=NULL) {
	lin <- tolower(gsub("^[[:blank:]]*|[[:blank:]]*$", "", readLines(file)))
	lin <- lin[nchar(lin) > 0]
	nexus <- any(grepl("charset", lin))
	iqtre <- any(grepl("modelfinder", lin))
	raxml <- any(grepl(",.+=.+-", lin))
	if (isTRUE(nexus)) {
		charsets <- sub("^charset[[:blank:]]+", "", grep("^charset", lin, value=TRUE))
		dfpart <- read_part_scheme(charsets)
		mb <- grepl("^partition", lin[max(grep("charset", lin)) + 1])
		iq <- any(grepl("charpartition", lin))
		if (isTRUE(mb)) {
			dfpart <- read_mb_model(lin, dfpart, param=param)	
		} else if (isTRUE(iq)) {
			dfpart <- read_iq_model(lin, dfpart, param=param)	
		} else {
			stop("The partition file was not recognized as either MrBayes or IQ-TREE nexus file.")
		}
	} else if (isTRUE(iqtre)) {
		dfpart <- read_iq_stdout(lin, param=param, name=name)
	} else if (isTRUE(raxml)) {
		dfpart <- read_part_scheme(lin)
	} else {
		stop("The partition file was not recognized as either nexus or RAxML-formatted partition file.")
	}
	return(dfpart)
}


read_part_scheme <- function(x) {
	getstring <- function(pattern, text) regmatches(text, regexpr(pattern, text))
	partlines <- grep(",.+=.+-|.+=.+-.+;", x)
	part <- x[partlines]
	part <- sub(";", "", part)
	part <- sub(" = ", "=", part)
	part <- gsub(",+", ",", gsub("[[:blank:]]", ",", part))
	type <- sub(",", "", getstring("^.+,", getstring("^.+=", part)))
	if (length(type) > 0) {
		part <- sub("^.+,", "", part)
		type <- sub("^DNA[[:alpha:]]", "DNA", type)
	}
	part <- strsplit(part, ",")
	npart <- length(part)
	if (npart == 0) {
		stop("There is no valid partition with nucleotide sequence data.")
	}
	part <- lapply(part, function(x) grep("^[DR]NA", x, invert=TRUE, value=TRUE))
	name <- unname(sapply(sapply(part, "[[", 1), function(x) sub("=", "", regmatches(x, regexpr("^.+=", x)))))
	name <- rep(name, sapply(part, length))
	part <- lapply(part, function(x) sub("^.+=", "", x))
	codon <- unlist(lapply(part, function(x) grepl("\\\\3", x)))
	part <- lapply(part, function(x) gsub("\\\\[[:digit:]]", "", x))
	start <- unlist(lapply(part, function(x) as.numeric(gsub("-.+$", "", x))))
	end <- unlist(lapply(part, function(x) as.numeric(gsub("^.+-", "", x))))	
	dfpart <- data.frame(name=name, start=start, end=end, codon=as.numeric(codon))
	codonpart <- unique(end[codon])
	for (i in seq_along(codonpart)) {
		ilocus <- dfpart$end == codonpart[i]
		dfpart$codon[ilocus] <- order(dfpart$start[ilocus])
		dfpart$start[ilocus] <- min(dfpart$start[ilocus])
	}
	dfpart$codon <- ifelse(dfpart$codon != 0, dfpart$codon, NA)
	return(dfpart)
}


read_iq_model <- function(x, dfpart, param) {
	getstring <- function(pattern, text) regmatches(text, gregexpr(pattern, text))
	equalfreqmodels <- c("jc", "jc69", "k80", "k2p", "tne", "k81", "k3p", "tpm2", "tpm3", "time", "tim2e", "tim3e", "tvme", "sym")
	npart <- sum(grepl("^charset", x))
	models <- gsub("[[:blank:]]", "", (x[grep("charpartition", x) + seq(npart)]))
	models <- sub(",$|;$", "", models)
	name <- gsub("^:|\\{.+\\}", "", getstring(":.+$", models))
	models <-  sub(":.+$", "", models)
	best_model <- grepl("\\{", x[grep("charpartition", x) + 1])
	if (isTRUE(best_model)) {
		if (isTRUE(param)) {
			values <- getstring("[[:alnum:]]+\\{[[:digit:]\\.,]+\\}", models)
		}
		models <- gsub("\\{[[:digit:]\\.,]+\\}", "", models)
	} else {
		param <- FALSE
	}
	gamma <- grepl("\\+g", models)
	free <- grepl("\\+r", models)
	rates <- ifelse(!gamma & !free, "Q", ifelse(free, "R", ifelse(gamma, "G", NA)))
	ncat <- sub("^[gr]", "", unlist(ifelse(gamma|free, getstring("[gr][[:digit:]]+", models[gamma|free]), NA)))
	pinv <- grepl("\\+i", models)
	eqfreq <- grepl("\\+fq", models)
	optfreq <- grepl("\\+fo", models)
	freq <- ifelse(optfreq, "FO", ifelse(eqfreq, "FQ", rep("F", length(models))))
	models <- sub("\\+.+$", "", models)
	freq <- ifelse(models %in% equalfreqmodels, "FQ", freq)
	dfmodels <- data.frame(model=toupper(models), rates=rates, ncat=ncat, pinv=pinv, freq=freq)
	if (isTRUE(param)) {
		dfmodels$subst <- sub("^[[:alnum:]]+", "", sapply(values, "[[", 1))
		values <- lapply(values, "[", -1)
		dfmodels$stfreq <- sub("^[[:alnum:]]+", "", sapply(values, function(v) unlist(getstring("f.*\\{[[:digit:]\\.,]+\\}", v))))
		dfmodels$gshape <- sub("^[[:alnum:]]+", "", sapply(values, function(v) unlist(getstring("g.*\\{[[:digit:]\\.,]+\\}", v))))
		dfmodels$gshape[dfmodels$gshape == "(0)"] <- NA
	}
	dfpart <- data.frame(dfpart, dfmodels[match(dfpart$name, name),], row.names=NULL)
	return(dfpart)
}


read_mb_model <- function(x, dfpart, param) {
	getstring <- function(pattern, text) regmatches(text, regexpr(pattern, text))
	x <- sub(";$", "", x)
	ord <- x[max(grep("charset", x)) + 1]
	ord <- sub("^.+:", "", gsub("[[:blank:]]", "", ord))
	ord <- unlist(strsplit(ord, ","))
	ord <- match(dfpart$name, ord)
	models <- grep("^lset.+nst", x, value=TRUE)
	which <- as.numeric(getstring("[[:digit:]]+", getstring("applyto=\\([[:digit:]]+\\)", models)))
	models <- sub("lset.+applyto=\\([[:digit:]]+\\).+nst", "nst", models)[which]
	model <- sub("nst=", "", getstring("nst=[126]", models))
	model <- unname(c("1"="JC", "2"="HKY", "6"="GTR")[model])
	gamma <- grepl("rates=gamma|rates=invgamma", models)
	rates <- ifelse(gamma, "G", "Q")
	pinv <- grepl("rates=propinv", models)
	ncat <- ifelse(gamma, 0, NA)
	ncat[gamma] <- as.numeric(getstring("[[:digit:]]+", getstring("ngammacat=[[:digit:]]+", models[gamma])))
	freq <- grep("prset.+statefreqpr=fixed", x, value=TRUE)
	if (length(freq) == 0) {
		freq <- rep("FO", length(models))
	} else if (grepl("applyto", freq)) {
		which <- as.numeric(getstring("[[:digit:]]+", getstring("applyto=\\([[:digit:]]+\\)", freq)))
		freq <- ifelse(grepl("equal", freq), "FQ", "F")
		freq <- ifelse(seq_along(models) %in% which, freq[order(which)], "FO")
	} else {
		freq <- rep(ifelse(grepl("equal", freq), "FQ", "F"), length(models))
	}
	dfmodels <- data.frame(model=model, rates=rates, ncat=ncat, pinv=pinv, freq=freq)
	dfpart <- data.frame(dfpart, dfmodels[ord,], row.names=NULL)
	return(dfpart)
}



read_iq_stdout <- function(x, param, name) {
	getstring <- function(pattern, text) regmatches(text, regexpr(pattern, text))
	equalfreqmodels <- c("jc", "jc69", "k80", "k2p", "tne", "k81", "k3p", "tpm2", "tpm3", "time", "tim2e", "tim3e", "tvme", "sym")
	if (any(grepl("partition", x))) {
		stop("The use of IQ-TREE standard output is restricted to unpartitioned data.")
	}
	name <- ifelse(is.null(name), NA, name)
	nsites <- as.numeric(getstring("^[[:digit:]]+", getstring("[[:digit:]]+ nucleotide site", x)))
	model <- sub("^.*:[[:blank:]]*", "", getstring("according to.+$", x))
	model <- unlist(strsplit(model, "\\+"))
	gamma <- any(grepl("^g", model))
	free <- any(grepl("^r", model))
	rates <- c("G", "R", "Q")[c(gamma, free, !gamma & !free)]
	ncat <- ifelse(gamma|free, getstring("[[:digit:]]+", model[gamma|free]), NA)
	pinv <- any( model== "i")
	eqfreq <- any(grepl("^fq", model)) | any(sapply(equalfreqmodels, grepl, x=model))
	optfreq <- any(grepl("^fo", model))
	freq <- c("FO", "FQ", "F")[c(optfreq, eqfreq, !optfreq & !eqfreq)]
	dfpart <- data.frame(name=name, start=1, end=nsites, codon=NA, model=toupper(model[1]), rates=rates, ncat=ncat, pinv=pinv, freq=freq)
	if (isTRUE(param)) {
		subst <- getstring("[[:digit:]\\.]+", sapply(c("a-c","a-g","a-t","c-g","c-t","g-t"), grep, x=x, value=TRUE))
		dfpart$subst <- paste0("{", paste(subst, collapse=","), "}")
		stfreq <- getstring("[[:digit:]\\.]+", sapply(c("pi\\(a\\)","pi\\(c\\)","pi\\(g\\)","pi\\(t\\)"), grep, x=x, value=TRUE))
		dfpart$stfreq <- paste0("{", paste(stfreq, collapse=","), "}")
		dfpart$gshape <- NA
	}
	return(dfpart)
}

