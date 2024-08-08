#' Principal component analysis.
#' 
#' @description
#' Principal component analysis of binary coded genotype data.
#' 
#' @param b matrix of binary coded genotype data, possibly the output of [make_binary_snps].
#' @returns A list with components `values` and `vectors` with results of the spectral decomposition of
#'   internally computed pairwise dissimilarity matrix. 
#' @details The method would be called 'principal coordinate analysis' or 'classical multidimensional scaling' by other authors
#'   as it is based on the decomposition of a dissimilarity matrix. The function is an implementation of approach described
#'   by Patterson et al. 2006 (PLoS Genet 2: e190). The spectral decomposition is performed by [base::eigen].
#' @export

pcomp <- function(b) {
	if (any(is.na(b))) {
		n <- nrow(b)
		X <- matrix(0, n, n)
		for (i in 1:n) {
			X[i,i] <- sum(b[i,]^2, na.rm=TRUE)
		}
		for (i in 1:(n-1)) {
			for (j in (i+1):n) {
				X[i,j] <- X[j,i] <- sum(b[i,] * b[j,], na.rm=TRUE)
			}
		}
		X <- X / n
	} else {
		X <- b %*% t(b) / nrow(b) 
	}
	eig <- base::eigen(X)
	dimnames(eig$vectors) <- list(rownames(b), paste0("PC", seq(nrow(X))))
	return(eig)	
}
