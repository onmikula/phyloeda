### FUNCTION
# seqlen
### DESCRIPTION
# calculates sequence lengths in the input alignment
### ARGUMENTS
# input: name of the file (fasta or phylip interleaved format) or sequence alignment in 'DNAbin' or 'matrix' format
### VALUE
# vector of sequence lengths
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
	return(rowSums(seq != "-"))
}
