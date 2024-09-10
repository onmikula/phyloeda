#' MrBayes input file.
#' 
#' @description
#' Writes input file for MrBayes analysis including data matrix and the block of MrBayes commands.
#'
#' @param seq a sequence alignment in `matrix` or `DNAbin` format or a name of the file with the alignment.
#' @param file a name of the MrBayes file input to be created by the function.
#' @param mcmc a vector of MCMC parameters, named by the parameter names used by MrBayes,
#'   the parameters not included are given default values (see Details).
#' @param sumt a vector of sumt command parameters, named by the parameter names used by MrBayes,
#'   the parameters not included are given default values (see Details). If `FALSE`, it is omitted.
#' @param prefix a stem of the name 
#' @param priors not implemented yet
#' @param topology 
#' @param interleaved logical, whether to write sequence data in interleaved format
#' @details The `mcmc` defaults are: ngen=2500000, samplefreq=1000, relburnin="yes", burninfrac=0.1, savebrlens="yes", nchains=4, nruns=4.
#'   The `sumt` defaults are: contype="allcompat", minpartfreq="0.0", calctreeprobs="no", showtreeprobs="no", hpd="yes", relburnin="yes",
#'   burninfrac=0.1.
#' @export

write_mrbayes_input <- function(seq, file, part, mcmc=NULL, sumt=NULL, prefix=NULL, priors=NULL, topology=NULL, interleaved=FALSE) {
	if (is.character(seq) & length(seq) == 1) {
		seq <- toupper(ape::read.dna(seq, format="fasta", as.matrix=TRUE, as.character=TRUE))
	}
	if (is.null(prefix)) {
		prefix <- gsub("^.*\\/|\\..*$", "", file)
	}
	mcmc <- c(mcmc, ngen=2500000, samplefreq=1000, relburnin="yes", burninfrac=0.1, savebrlens="yes", nchains=4, nruns=4, filename=prefix)
	names(mcmc) <- tolower(names(mcmc))
	mcmc <- mcmc[!duplicated(names(mcmc))]
	if (isFALSE(sumt)) {
		sumt <- NULL
	} else {
		sumt <- c(sumt, contype="allcompat", minpartfreq="0.0", calctreeprobs="no", showtreeprobs="no", hpd="yes", relburnin="yes", burninfrac=0.1)
		names(sumt) <- tolower(names(sumt))
		sumt <- sumt[!duplicated(names(sumt))]
	}
	if (missing(part)) {
		if (missing(seq)) {
			part <- data.frame(name=NA, start=NA, end=NA, codon=NA, model="HKY", rates="gamma", ncat=4, pinv=FALSE, freq=FALSE)
		} else {
			part <- data.frame(name=NA, start=1, end=ncol(seq), codon=NA, model="HKY", rates="gamma", ncat=4, pinv=FALSE, freq=FALSE)		
		}
	} else if (is.character(part) & length(part) == 1) {
		part <- read_part(part)
	}
	if (any(part$rates == "R")) {
		part$rates[part$rates == "R"] <- "G"
		warning("Free rate variation model (R) was changed to gamma (G).")
	}
	mbblock <- prepare_mbblock(part=part, mcmc=mcmc, sumt=sumt, priors=priors, topology=topology)
#	mbinput <- c(capture.output(ape::write.nexus.data(seq, file=stdout(), interleaved=TRUE, charsperline=80, gap="!", missing="-")), rep("\n", 2), mbblock)
	ape::write.nexus.data(seq, file=file, interleaved=interleaved, charsperline=80, gap="-", missing="N")
	mbinput <- c(readLines(file), "\n", mbblock)
	writeLines(mbinput, file)
}



prepare_mbblock <- function(part, mcmc, sumt, priors=NULL, topology=NULL) {
	
	if (nrow(part) > 1 & any(is.na(part$name))) {
		stop("If multiple partitions are defined, all of them must be named.")
	}
	npart <- length(unique(part$name))
	if (npart == 1) {
		charsets <- NULL
		partitions <- NULL
	} else {
		charsets <- part[,c("name", "start", "end", "codon")]
		bycodon <- all(!is.na(charsets$codon)) & all(charsets$start == charsets$start[1])
		charsets$start <- ifelse(!is.na(charsets$codon), charsets$start + charsets$codon - 1, charsets$start)
		charsets$codon <- ifelse(!is.na(charsets$codon), "\\3", "")
		charsets$range <- paste0(charsets$start, "-", charsets$end, charsets$codon)
		for (i in which(duplicated(charsets$name))) {
			ii <- which(charsets$name == charsets$name[i])
			charsets$range[ii[1]] <- paste(charsets$range[ii], collapse=" ")
		}
		charsets <- charsets[!duplicated(charsets$name),]
		partitions <- paste("partition", c("part","bycodon")[bycodon+1], "=", nrow(charsets))
		partitions <- paste0(partitions, ": ", paste(charsets$name, collapse=", "), ";")
		partitions <- c(partitions, paste0("set partition=", c("part","bycodon")[bycodon+1], ";"))
		charsets <- paste0(paste("charset", charsets$name, "=", charsets$range), ";")
		part <- part[!duplicated(part$name),]
		npart <- nrow(part)
	}

	nst <- c(JC=1, JC69=1, F81=1, K80=2, K2P=2, HKY=2, HKY85=2, SYM=6, GTR=6)[part$model]
	if (npart == 1) {
		lset <- list("lset ")
	} else {
		lset <- as.list(paste0("lset applyto=(", seq(npart), ") "))
	}
	for (i in seq(npart)) {
		lset[[i]] <- paste0(lset[[i]], "nst=", nst[i])
		if (isTRUE(part$pinv[i]) & part$rates[i] == "Q") {
			lset[[i]] <- paste(lset[[i]], "rates=propinv;")
		} else if (isTRUE(part$pinv[i]) & part$rates[i] == "G") {
			lset[[i]] <- paste(lset[[i]], "rates=invgamma", paste0("ngammacat=", part$ncat[i], ";"))
		} else if (isFALSE(part$pinv[i]) & part$rates[i] == "G") {
			lset[[i]] <- paste(lset[[i]], "rates=gamma", paste0("ngammacat=", part$ncat[i], ";"))
		} else {
			lset[[i]] <- paste0(lset[[i]], ";")
		}
	}
	lset <- unlist(lset)

	unlink <- NULL
	if (npart > 1) {
		if (sum(part$freq == "FO") > 1) unlink <- c(unlink, paste0("statefreq=(", paste(which(part$freq == "FO"), collapse=","), ")"))
		if (sum(nst == 2) > 1) unlink <- c(unlink, paste0("tratio=(", paste(which(nst == 2), collapse=","), ")"))
		if (sum(nst == 6) > 1) unlink <- c(unlink, paste0("revmat=(", paste(which(nst == 6), collapse=","), ")"))
		if (sum(part$pinv) > 1) unlink <- c(unlink, paste0("pinvar=(", paste(which(part$pinv), collapse=","), ")"))
		if (sum(part$rates != "Q") > 1) unlink <- c(unlink, paste0("shape=(", paste(which(part$rates != "Q"), collapse=","), ")"))
		unlink <- paste0(paste(c("unlink", unlink), collapse=" "), ";")
	}

	prset <- NULL
	if (any(part$freq == "FQ"))
		prset <- c(prset, paste0("prset applyto=(", paste(which(part$freq == "FQ"), collapse=","), ") statefreqpr=fixed(equal);"))
	if (any(part$freq == "F")) {
		stfreq <- gsub("\\{|\\}", "", part$stfreq[part$freq == "F"])
		prset <- c(prset, paste0("prset applyto=(", which(part$freq == "F"), ") statefreqpr=fixed(", stfreq, ");"))
	}
	if (npart > 1) prset <- c(prset, "prset applyto=(all) ratepr=variable;")
	if (!is.null(topology)) prset <- c(prset, "prset topologypr=fixed(topology);")

	mcmcp <- paste0(paste(c("mcmcp", paste(names(mcmc), mcmc, sep="=")), collapse=" "), ";")
	if (!is.null(sumt)) {
		sumt <- paste0(paste(c("sumt", paste(names(sumt), sumt, sep="=")), collapse=" "), ";")
	}

	if (!is.null(topology)) {
		topology <- read_tree(topology)
		topology$edge.length <- NULL
		topology$node.label <- NULL
		topology <- ape::write.tree(topology)
		trees <- c("begin trees;", paste("tree", "topology", "=", topology), "end;", "")
	} else {
		trees <- NULL
	}

	mrbayes <- c(trees, "begin mrbayes;", charsets, partitions, lset, unlink, prset, priors, mcmcp, "mcmc;", sumt, "end;")
	
	return(mrbayes)

}

