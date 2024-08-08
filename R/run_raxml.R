#' RAxML analysis.
#' 
#' @description
#' Runs RAxML analysis specified by the arguments, including pre- and post-processing of the files.
#' 
#' @param command a character string with the system command. If not supplied, either "raxmlHPC-SSE3" or "raxmlHPC.exe" is used,
#'   depending on the platform. Also the path to directory with the appropriate binary file can be specified
#'   instead of the command itself.
#' @param args a character vector in the form accepted by `args` argument of [base::system2] or a data frame with arguments
#'   to `command` (argument names in the 1st column, values in the 2nd column) or a name of a tab-delimited file with the arguments.
#' @param prefix a character string to be included in names of all output files. It has priority over the specification
#'   in `args` argument. If not specified, whether through `prefix` or `args`, the input file name without extension is used.
#' @param outdir a character string with the path to output directory. It has priority over the path included in `prefix`.
#'   If not specified in any way, it is assumed to be the same as the working directory.
#' @param retain regular expressions with identifiers of outputs to be retained. If specified, the other outputs are deleted.
#' @param stdout,stderr where output to ‘stdout’ or ‘stderr’ should be sent, as in [base::system2].
#' @param seed random seed for the analysis ("-p" argument of RAxML). It has priority over the specification in `args`.
#' @details The minimum `args` contains "-s" and the name of the input file. The subsotution model ("-m") is "GTRGAMMA" by default. 
#'   The input directory is either retrieved from the input name (if it contains a path) or assummed to be the working directory.
#'   The argument names must be supplied  with the leading dash ("-"). If `args` is a data frame, the arguments with double dash ("--")
#'   have `NA` in the 2nd column.
#'   For the reference to RAxML commands, consult https://cme.h-its.org/exelixis/web/software/raxml/hands_on.html
#'   or RAxML manual at https://cme.h-its.org/exelixis/resource/download/NewManual.pdf.
#' @export

run_raxml <- function(command, args, prefix, outdir, retain, stdout="", stderr="", seed=123) {

	if (missing(command)) {
		command <- ifelse(.Platform$OS.type == "unix", "raxmlHPC-SSE3", "raxmlHPC.exe")
	}
	if (dir.exists(command)) {
		command <- paste0(sub("\\/$", "", command), ifelse(.Platform$OS.type == "unix", "/raxmlHPC-SSE3", "/raxmlHPC.exe"))
	}

	if (missing(args)) {
		stop("The 'args' argument was not specified with no default.")
	}
	if (is.character(args) & length(args) == 1) {
		args <- read.delim(args)
	}
	if (!any(as.character(args) == "-s")) {
		stop("The sequence file is not specified in 'args'.")
	}
	istable <- is.data.frame(args) | is.matrix(args)
	dupl <- ifelse(istable, any(duplicated(args[,1])), any(duplicated(args[seq_along(args)%%2 == 1])))
	if (isTRUE(dupl)) {
		stop("'args' contains duplicated entries.")
	}
	if (istable) {
		args <- na.omit(unname(unlist(Map(c, as.character(args[,1]), as.character(args[,2])))))
	}

	if (missing(prefix)) {
		if (any(args == "-n")) {
			prefix <- args[which(args == "-n") + 1]
		} else {
			prefix <- sub("\\.[[:alnum:]]*$", "", args[which(args == "-s") + 1])
		}
	}
	if (missing(outdir)) {
		if (grepl(".+\\/", prefix)) {
			outdir <- regmatches(prefix, regexpr("^.+\\/", prefix))
		} else {
			outdir <- getwd()
		}
	}
	prefix <- sub("\\/", "", prefix)
	if (any(args == "-n")) {
		args[which(args == "-n") + 1] <- prefix
	} else {
		args <- c(args, "-n", prefix)
	}
	if (!any(args == "-m")) {
		args <- c(args, "-m", "GTRGAMMA")
	}
	if (any(args == "-p")) {
		args[which(args == "-p") + 1] <- seed
	} else {
		args <- c(args, "-p", seed)
	}

			
	seqfile <- args[which(args == "-s") + 1]
	partfile <- ifelse(any(args == "-q"), args[which(args == "-q") + 1], NA)
	treefile <- ifelse(any(args == "-t"), args[which(args == "-t") + 1], NA)
	seqout <- grepl("\\/", seqfile)
	partout <- grepl("\\/", partfile)
	treeout <- grepl("\\/", treefile)
	if (seqout) {
		cseqfile <- sub("^.*\\/", "", seqfile)
		file.copy(seqfile, cseqfile)
		args[which(args == "-s") + 1] <- cseqfile
	} else {
		cseqfile <- NULL
	}
	if (partout) {
		cpartfile <- sub("^.*\\/", "", partfile)
		file.copy(partfile, cpartfile)
		args[which(args == "-q") + 1] <- cpartfile
	} else {
		cpartfile <- NULL
	}
	if (treeout) {
		ctreefile <- sub("^.*\\/", "", treefile)
		file.copy(treefile, ctreefile)
		args[which(args == "-t") + 1] <- ctreefile
	} else {
		ctreefile <- NULL
	}

	system2(command, args=args, stdout=stdout, stderr=stderr)

	cfiles <- unlist(list(cseqfile, cpartfile, ctreefile)[c(seqout, partout, treeout)])
	if (length(cfiles) > 0) {
		invisible(file.remove(cfiles))
	}
	invisible(file.remove(paste0(args[which(args == "-s") + 1], ".reduced")))
	outputs <- list.files(".")
	fixed <- any(sapply(unlist(strsplit(prefix, "")), grepl, x=".|()[{^$*+?⁠", fixed=TRUE))
	outputs <- grep(prefix, outputs, value=TRUE, fixed=fixed)
	outputs <- grep("^RAxML_", outputs, value=TRUE, fixed=FALSE)
	outputs <- setdiff(outputs, c(seqfile, partfile, treefile))
	if (missing(retain)) {
		retained <- outputs
	} else {
		retained <- unlist(lapply(retain, grep, outputs, value=TRUE))
		retained <- sort(unique(c(retained, unlist(lapply(retain, grep, outputs, value=TRUE, fixed=TRUE)))))
	}
	if (outdir != getwd()) {
		invisible(file.copy(retained, paste0(sub("\\/+$","", outdir), "/", retained)))
		invisible(file.remove(outputs))
	} else {
		invisible(file.remove(setdiff(outputs, retained)))
	}

}



#' Evolutionary placement analysis.
#' 
#' @description
#' Runs evolutionary placement analysis (EPA) under RAxML, equivalent to `run_raxml` with appropriately set arguments.
#' 
#' @param command a character string with the system command. If not supplied, either "raxmlHPC-SSE3" or "raxmlHPC.exe" is used,
#'   depending on the platform. Also the path to directory with the appropriate binary file can be specified
#'   instead of the command itself.
#' @param args as in [run_raxml], but not only "-s" but also "-t" specification is obligatory and "-f" is fixed to "v".
#'   For defaults of other parameters, check `Details`.
#' @param prefix a character string to be included in names of all output files. It has priority over the specification
#'   in `args` argument. If not specified, whether through `prefix` or `args`, the input file name without extension is used.
#' @param outdir a character string with the path to output directory. It has priority over the path included in `prefix`.
#'   If not specified in any way, it is assumed to be the same as the working directory.
#' @param retain regular expressions with identifiers of outputs to be retained. The .jplace output is always retained, but
#'   others can be retained as well (specify `retain="^RAxML_"` to retain all fo them).
#'   changed, the other outputs than the ".jplace" file are deleted.
#' @param rename if `TRUE` (default) the name of .jplace output file is changed to `paste0(prefix, ".jplace")`. 
#'   Alternatively, the name can be specified by a character string or left as it is (`rename=FALSE`).
#' @param stdout,stderr where output to ‘stdout’ or ‘stderr’ should be sent, as in [base::system2].
#' @details The minimum `args` contains "-s" and the name of the input file as well as "-t" and the name of the tree file.
#'   The input directory is either retrieved from the input name (if it contains a path) or assummed to be the working directory.
#'   Internally, the function always adds `c("-f", "v")` to `args` and if not specified otherwise by the user, also adds defaults
#'   for EPA parameters (`"--epa-keep-placements=100"` and `"--epa-prob-threshold=0.01"`) and for convenience also "--HKY85"
#'   as default substitution model. For the reference to RAxML commands, consult RAxML manual at
#'   https://cme.h-its.org/exelixis/resource/download/NewManual.pdf.
#' @export

run_epa <- function(command, args, prefix, outdir, retain="\\.jplace$", rename=TRUE, stdout="", stderr="") {

	if (missing(command)) {
		command <- ifelse(.Platform$OS.type == "unix", "raxmlHPC-SSE3", "raxmlHPC.exe")
	}
	if (dir.exists(command)) {
		command <- paste0(sub("\\/$", "", command), ifelse(.Platform$OS.type == "unix", "/raxmlHPC-SSE3", "/raxmlHPC.exe"))
	}

	if (missing(args)) {
		stop("The 'args' argument was not specified with no default.")
	}
	if (is.character(args) & length(args) == 1) {
		args <- read.delim(args)
	}
	if (!any(as.character(args) == "-s")) {
		stop("The sequence file is not specified in 'args'.")
	}
	if (!any(as.character(args) == "-t")) {
		stop("The tree file is not specified in 'args'.")
	}
	istable <- is.data.frame(args) | is.matrix(args)
	dupl <- ifelse(istable, any(duplicated(args[,1])), any(duplicated(args[seq_along(args)%%2 == 1])))
	if (isTRUE(dupl)) {
		stop("'args' contains duplicated entries.")
	}
	if (istable) {
		args <- na.omit(unname(unlist(Map(c, as.character(args[,1]), as.character(args[,2])))))
	}
	if (missing(prefix)) {
		if (any(args == "-n")) {
			prefix <- args[which(args == "-n") + 1]
		} else {
			prefix <- sub("\\.[[:alnum:]]*$", "", args[which(args == "-s") + 1])
		}
	}
	if (missing(outdir)) {
		if (grepl(".+\\/", prefix)) {
			outdir <- regmatches(prefix, regexpr("^.+\\/", prefix))
		} else {
			outdir <- getwd()
		}
	}
	prefix <- sub("\\/", "", prefix)
	if (any(args == "-n")) {
		args[which(args == "-n") + 1] <- prefix
	} else {
		args <- c(args, "-n", prefix)
	}
	if (!any(args == "-m")) {
		args <- c(args, "-m", "GTRGAMMA")
	}
	if (any(args == "-f")) {
		args[which(args == "-f") + 1] <- "v"
	} else {
		args <- c("-f", "v", args)
	}
	epa <- c(keep="--epa-keep-placements", prob="--epa-prob-threshold", cum="--epa-accumulated-threshold")
	hasepa <- sapply(lapply(epa, grepl, x=args), any)
	if (any(hasepa[1:2])) {
		if (hasepa[3]) {
			args <- args[-grep("--epa-accumulated-threshold", args)]
		}
		hasepa[3] <- TRUE
	} else if (!any(hasepa)) {
		hasepa[3] <- TRUE
	}
	args <- c(args, c("--epa-keep-placements=100", "--epa-prob-threshold=0.01", "--epa-accumulated-threshold=0.99")[!hasepa])
	if (!any(c("--JC69", "--K80") %in% args)) {
		args <- c(args, "--HKY85")
	}
			
	seqfile <- args[which(args == "-s") + 1]
	partfile <- ifelse(any(args == "-q"), args[which(args == "-q") + 1], NA)
	treefile <- ifelse(any(args == "-t"), args[which(args == "-t") + 1], NA)
	seqout <- grepl("\\/", seqfile)
	partout <- grepl("\\/", partfile)
	treeout <- grepl("\\/", treefile)
	if (seqout) {
		cseqfile <- sub("^.*\\/", "", seqfile)
		file.copy(seqfile, cseqfile)
		args[which(args == "-s") + 1] <- cseqfile
	} else {
		cseqfile <- NULL
	}
	if (partout) {
		cpartfile <- sub("^.*\\/", "", partfile)
		file.copy(partfile, cpartfile)
		args[which(args == "-q") + 1] <- cpartfile
	} else {
		cpartfile <- NULL
	}
	if (treeout) {
		ctreefile <- sub("^.*\\/", "", treefile)
		file.copy(treefile, ctreefile)
		args[which(args == "-t") + 1] <- ctreefile
	} else {
		ctreefile <- NULL
	}

	system2(command, args=args, stdout=stdout, stderr=stderr)

	cfiles <- unlist(list(cseqfile, cpartfile, ctreefile)[c(seqout, partout, treeout)])
	if (length(cfiles) > 0) {
		invisible(file.remove(cfiles))
	}
	invisible(file.remove(paste0(args[which(args == "-s") + 1], ".reduced")))
	outputs <- list.files(".")
	fixed <- any(sapply(unlist(strsplit(prefix, "")), grepl, x=".|()[{^$*+?⁠", fixed=TRUE))
	outputs <- grep(prefix, outputs, value=TRUE, fixed=fixed)
	outputs <- grep("^RAxML_", outputs, value=TRUE, fixed=FALSE)
	outputs <- setdiff(outputs, c(seqfile, partfile, treefile))
	
	retain <- unique(c("\\.jplace$", retain))
	retained <- unlist(lapply(retain, grep, outputs, value=TRUE))
	retained <- sort(unique(c(retained, unlist(lapply(retain, grep, outputs, value=TRUE, fixed=TRUE)))))
	if (isTRUE(rename)) {
		jplacefile <- paste0(prefix, ".jplace")
	} else if (is.character(rename)) {
		jplacefile <- rename
	}
	if (!isFALSE(rename)) {
		copies <- ifelse(grepl("\\.jplace$", retained), jplacefile, retained)
	} else {
		copies <- retained
	}
	if (outdir != getwd()) {	
		copies <- paste0(sub("\\/+$","", outdir), "/", copies)
	}
	renamed <- retained != copies
	if (any(renamed)) {
		invisible(file.rename(retained[renamed], copies[renamed]))		
	}
	removed <- setdiff(outputs, retained)
	if (length(removed) > 0) {
		invisible(file.remove(removed))	
	}

}
