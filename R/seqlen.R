#' Sequence lengths.
#' 
#' @description
#' Calculates lengths of sequences in the input alignment.
#'
#' @param input a character string, the name of file in .fasta or .phylip format
#'   or a sequence alignment in `DNAbin` or `matrix` format.
#' @returns An numeric vector with sequence lengths.
#' @export

seqlen <- function(input) {
	if (inherits(input, "DNAbin")) {
		seq <- toupper(as.character(input))
	} else if (inherits(input, "matrix")) {
		seq <- toupper(input)
	} else {
		if (substr(readLines(input, n=1), 1, 1) == ">") {
			seq <- toupper(ape::read.dna(input, format="fasta", as.character=TRUE, as.matrix=TRUE))
		} else {
			seq <- toupper(ape::read.dna(input, as.character=TRUE, as.matrix=TRUE))
		}	
	}
	seq <- ifelse(seq == "-" | seq == "N" | seq == "n" | seq == "?", "-", seq)
	len <- rowSums(seq != "-")
	return(len)
}
