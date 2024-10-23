#' Principal coordinate analysis.
#' 
#' @description
#' Principal coordinate analysis, a.k.a. classical multidimensional scaling.
#' 
#' @param d matrix of pairwise dissimilarities, possibly calculated from binary coded genotype data by [dissim.snps].
#' @returns A list with components `values` and `vectors` with results of the spectral decomposition of
#'   internally computed pairwise inter-individual covariance matrix. 
#' @details The function is an implementation of PCoA as described, e.g., in Numerical Ecology (Legendre et al. ----).
#'   The spectral decomposition is performed by [base::eigen].
#' @export

pcoord <- function(dst) {
	if (inherits(dst, "dist")) {
		d <- as.matrix(d)
	}
	n <- nrow(d)
	C <- diag(n) - (1 / n) * rep(1, n) %*% t(rep(1, n))
	A <- -0.5 * d^2
	delta <- C %*% A %*% C
	eig <- eigen(delta)
	pos <- eig$values >= 0
	val <- eig$values[pos]
	vec <- eig$vectors[,pos] %*% diag(sqrt(val))
	dimnames(vec) <- list(rownames(d), paste0("PCo", seq(sum(pos))))
	result <- list(vectors=vec, values=val)
	class(result) <- "pcoa"
	return(result)
}
