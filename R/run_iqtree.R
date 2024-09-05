#' IQ-TREE analysis.
#' 
#' @description
#' Runs IQ-TREE analysis specified by the arguments, including post-processing of the output files.
#' 
#' @param command a character string with the system command. If not supplied, either "iqtree2" or "iqtree.exe" is used,
#'   depending on the platform. Also the path to directory with the appropriate binary file can be specified
#'   instead of the command itself.
#' @param args a character vector in the form accepted by `args` argument of [base::system2] or a data frame with arguments
#'   to `command` (argument names in the 1st column, values in the 2nd column) or a name of a tab-delimited file with the arguments.
#' @param prefix a character string with prefix for all output files. It has priority over the specification
#'   through `args` argument. If not specified, whether through `prefix` or `args`, the input file name
#'   without extension is used. 
#' @param outdir a character string with the path to output directory. It has priority over the path included in prefix.
#'   If not specified in any way, it is assumed to be the same as the working directory.
#' @param retain regular expressions with identifiers of outputs to be retained. If specified, the other outputs are deleted.
#' @param coding logical, if `TRUE` a proper partition file is created and used in the substitution model selection.
#' @param stdout,stderr where output to ‘stdout’ or ‘stderr’ should be sent, as in [base::system2].
#' @details The minimum `args` contains single row with "-s" and the name of the input file. The input directory is
#'   either retrieved fron the input name (if it contains a path) or assummed to be the working directory. The argument
#'   names must be supplied with the leading dashes ("-" or "--", what is appropriate). For the reference to IQ-TREE commands,
#'   consult ... http://www.iqtree.org/doc/Command-Reference ...
#' @export

run_iqtree <- function(command, args, prefix, outdir, retain, coding=FALSE, stdout="", stderr="") {

	if (missing(command)) {
		command <- ifelse(.Platform$OS.type == "unix", "iqtree2", "iqtree.exe")
	}
	if (dir.exists(command)) {
		command <- paste0(sub("\\/$", "", command), ifelse(.Platform$OS.type == "unix", "/iqtree2", "/iqtree.exe"))
	}

	if (missing(args)) {
		stop("The 'args' argument was not specified with no default.")
	}
	if (is.character(args) & length(args) == 1) {
		args <- read.delim(args)
	}
	sarg <- ifelse(is.data.frame(args), any(args[,1] == "-s"), any(as.character(args) == "-s"))
	if (!sarg) {
		stop("The sequence file is not specified in 'args'.")
	}
	istable <- is.data.frame(args) | is.matrix(args)
	dupl <- ifelse(istable, any(duplicated(args[,1])), any(duplicated(args[seq_along(args)%%2 == 1])))
	if (isTRUE(dupl)) {
		stop("'args' contains duplicated entries.")
	} else if (istable) {
		args <- unname(unlist(Map(c, as.character(args[,1]), as.character(args[,2]))))
	}
	coding <- isTRUE(coding) & all(args != "-p")
	if (isTRUE(coding)) {
		partfile <- "temp_part.txt"
		nbp <- gsub("^\\s+|\\s+$", "", readLines(args[which(args == "-s") + 1], n=10))
		nbp <- nbp[nchar(nbp) > 0]
		if (grepl("^[[:digit:]]+\\s+[[:digit:]]+$", nbp[1])) {
			nbp <- as.numeric(regmatches(nbp[1], regexpr("[[:digit:]]+$", nbp[1])))
		} else {
			nbp <- ncol(read.fasta(args[which(args == "-s") + 1]))
		}
		args <- c(args, "-p", partfile)
		if (any(args == "-m")) {
			m <- which(args == "-m") + 1
			if (!grepl("MERGE", args[m])) {
				if (grepl("TESTNEW", args[m])) {
					args[m] <- sub("TESTNEW", "TESTNEWMERGE", args[m])
				} else if (grepl("TEST", args[m])) {
					args[m] <- sub("TEST", "TESTMERGE", args[m])
				} else if (grepl("MFP", args[m])) {
					args[m] <- sub("MFP", "MFP+MERGE", args[m])
				}
			}
		} else {
			args <- c(args, "-m", "TESTMERGE")
		}
		write_codon_part(scheme="1-2-3", file=partfile, first=1, last=nbp, name="pos", format="raxml")
	}
	if (missing(prefix)) {
		if (any(args == "--prefix")) {
			prefix <- args[which(args == "--prefix") + 1]
		} else {
			prefix <- sub("\\.[[:alnum:]]*$", "", args[which(args == "-s") + 1])
		}
	}
	if (missing(outdir)) {
		if (grepl(".+\\/", prefix)) {
			outdir <- regmatches(prefix, regexpr("^.+\\/", prefix))
		} else if (grepl(".+\\/", args[which(args == "-s") + 1])) {
			input <- args[which(args == "-s") + 1]
			outdir <- regmatches(input, regexpr("^.+\\/", input))	
		} else {
			outdir <- getwd()
		}
	}
	outdir <- ifelse(!grepl("\\/$", outdir), paste0(outdir, "/"), outdir)
	prefix <- paste0(outdir, sub("^.*\\/", "", prefix))
	if (any(args == "--prefix")) {
		args[which(args == "--prefix") + 1] <- prefix
	} else {
		args <- c(args, "--prefix", prefix)
	}
		
	system2(command, args=args, stdout=stdout, stderr=stderr)
	
	if (!missing(retain)) {
		outputs <- list.files(outdir, full.names=TRUE)
		prefix <- sub("^.*\\/", "", prefix)
		fixed <- any(sapply(unlist(strsplit(prefix, "")), grepl, x=".|()[{^$*+?⁠", fixed=TRUE))
		outputs <- grep(prefix, outputs, value=TRUE, fixed=fixed)
		retained <- unlist(lapply(retain, grep, outputs, value=TRUE))
		retained <- sort(unique(c(retained, unlist(lapply(retain, grep, outputs, value=TRUE, fixed=TRUE)))))
		retained <- c(args[which(args == "-s") + 1], retained)
		invisible(file.remove(setdiff(outputs, retained)))
	}	
	if (isTRUE(coding)) {
		invisible(file.remove(partfile))
	}

}

