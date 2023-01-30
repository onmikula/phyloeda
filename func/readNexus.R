# Function:
#	readNexus: reads sequence from (partitioned) .nexus file
# Arguments:
#	nex: name of the partitioned .nexus file
#	toupper: whether to export sequences in upper case letters
#	as.DNAbin: whether to export sequences as 'DNAbin' objects
# Value:
#	(a list of) alignment(s) in matrix or DNAbin format (sigle alignment if no partitions are defined in the .nexus file)
# Depends on: ape

readNexus <- function(nex, toupper=TRUE, as.DNAbin=FALSE) {
	lin <- readLines(nex)
	lin <- lin[nchar(lin) > 0]
	start <- which(lin == "MATRIX") + 1
	end <- base::intersect(which(lin == ";"), start:length(lin))[1] - 1
	dat <- unlist(strsplit(lin[start:end], "\\s+"))
	seq <- do.call(rbind, strsplit(dat[seq(length(dat)) %% 2 == 0], ""))
	rownames(seq) <- dat[seq(length(dat)) %% 2 == 1]
	last <- lin[(end + 1):length(lin)]
	part <- last[grepl("charset", last, ignore.case=TRUE)]
	if (length(part) > 0) {
		part <- gsub("charset|[[:blank:]]|;", "", part, ignore.case=TRUE)
		part <- strsplit(part, "=")
		ind <- lapply(strsplit(sapply(part, "[", 2), "-"), as.numeric)
		seq <- lapply(by(t(seq), rep(seq_along(ind), sapply(ind, diff) + 1), identity), t)
		names(seq) <- sapply(part, "[", 1)
		if (isTRUE(toupper) & isFALSE(as.DNAbin)) {
			seq <- lapply(seq, toupper)
		}
		if (isTRUE(as.DNAbin)) {
			seq <- lapply(seq, ape::as.DNAbin)
		}
	} else {
		if (isTRUE(toupper) & isFALSE(as.DNAbin)) {
			seq <- toupper(seq)
		}
		if (isTRUE(as.DNAbin)) {
			seq <- ape::as.DNAbin(seq)
		}
	}
	return(seq)
}


# Function:
#	readNexusMatrix: reads sequence from unpartitioned .nexus file with sequences in matrix format
# Arguments:
#	nex: name of the .nexus file

readNexusMatrix <- function(nex) {
	x <- readLines(nex)
	x <- gsub("^[[:blank:]]*", "", x)
	x <- x[nchar(x) > 0]
	start <- grep("^matrix", x) + 1
	end <- grep("end;$", x, ignore.case=TRUE) - 2
	mat <- do.call(rbind, strsplit(x[start:end], "\\s"))
	if (any(grepl("interleave", x))) {
		mat <- lapply(split(mat[,2], mat[,1]), paste, collapse="")
		seqdata <- do.call(rbind, strsplit(unlist(mat), ""))
		rownames(seqdata) <- names(mat)	
	} else {
		seqdata <- do.call(rbind, strsplit(mat[,2], ""))
		rownames(seqdata) <- mat[,1]	
	}
	return(seqdata)
}
