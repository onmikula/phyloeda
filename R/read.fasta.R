#' Reading fasta file.
#' 
#' @description
#' A wrapper for [ape::read.dna] with options `format="fasta"`, `as.matrix=TRUE` and `as.character=TRUE`.
#'
#' @param file character string, the name of the input fasta file.
#' @returns A matrix of mode `character` with nucleotides represented by upper case letters.
#'   The coding of gap, missing or ambiguous data is left as it was in the input file.
#' @export

read.fasta <- function(file) toupper(ape::read.dna(file, format="fasta", as.matrix=TRUE, as.character=TRUE))



#' Writing fasta file.
#' 
#' @description
#' A wrapper for [ape::write.dna] with options `format="fasta"`, `as.matrix=TRUE` and `as.character=TRUE`.
#'
#' @param seq a matrix of mode `character` or an object of class `DNAbin`.
#' @param file character string, the ame of the output fasta file.
#' @details It writes nucleotide in upper case, irrespective of their original coding.
#'   The coding of gap, missing or ambiguous data is left as it was in the input object.
#' @export

write.fasta <- function(seq, file) ape::write.dna(toupper(as.matrix(seq)), file, format="fasta")
