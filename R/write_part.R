#' Write partition file.
#' 
#' @description
#' Writes nexus or RAxML-formatted partition files.
#'
#' @param part a data frame with information about partition scheme, formatted as an output of [read_part].
#' @param file a name of the partition file, if missing, the output is printed to console.
#' @param format a character string, either `"raxml"` or `"nexus"`.
#' @param type currently, only `"DNA"` is assumed here.
#' @export

write_part <- function(part, file, format=c("raxml","nexus"), type="DNA") {
	if (missing(file)) {
		file <- stdout()
	}
	format <- format[1]
	charsets <- part[,c("name", "start", "end", "codon")]
	bycodon <- all(!is.na(charsets$codon)) & all(charsets$start == charsets$start[1])
	charsets$start <- ifelse(!is.na(charsets$codon), charsets$start + charsets$codon - 1, charsets$start)
	charsets$codon <- ifelse(!is.na(charsets$codon), "\\3", "")
	charsets$range <- paste0(charsets$start, "-", charsets$end, charsets$codon)
	for (i in which(duplicated(charsets$name))) {
		ii <- which(charsets$name == charsets$name[i])
		charsets$range[ii[1]] <- paste(charsets$range[ii], collapse=c(raxml=", ", nexus=" ")[format])
	}
	charsets <- charsets[!duplicated(charsets$name),]
	if (format == "nexus") {
		charsets <- paste0(paste("charset", charsets$name, "=", charsets$range), ";")
		charsets <- c("#nexus", "begin sets;", charsets, "end;")
	} else if (format == "raxml") {
		charsets <- paste(paste0(type, ","), charsets$name, "=", charsets$range)
	}
	writeLines(charsets, file)
}




#' Codon partition scheme.
#' 
#' @description
#' Writes nexus or RAxML-formatted partition file for a single codon-partitioned locus.
#'
#' @param scheme a character string symbolically representing partitioning scheme.
#' @param file a name of the partition file, if missing, the output is printed to console.
#' @param last numeric, number of sites in the alignement.
#' @param first numeric, position of the first nucleotide (within a longer genomic sequence).
#' @param name a character string, label(s) used to name partitions, if necessary, it is recycled and given the codon position numbers.
#' @param format a character string, either `"raxml"` or `"nexus"`.
#' @export

write_codon_part <- function(scheme, file, last, first=1, name="pos", format=c("raxml","nexus")) {
	if (missing(file)) {
		file <- stdout()
	}
	format <- format[1]
	codon <- setNames(paste0(first - 1 + 1:3, "-", first - 1 + last, "\\3"), paste0(name, 1:3))
	scheme <- lapply(sapply(unlist(strsplit(scheme, "-")), strsplit, split=""), as.numeric)
	charsets <- lapply(scheme, function(ii) codon[ii])
	charsets <- sapply(charsets, function(x) paste(paste(names(x), collapse="_"), "=", paste(x, collapse=c(raxml=", ", nexus=" ")[format])))
	if (format == "nexus") {
		charsets <- paste0(paste("charset", charsets), ";")
		charsets <- c("#nexus", "begin sets;", charsets, "end;")
	} else if (format == "raxml") {
		charsets <- paste("DNA,", charsets)
	}
	writeLines(charsets, file)
}
