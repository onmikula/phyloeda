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
	eig <- eigen(X)
	dimnames(eig$vectors) <- list(rownames(b), paste0("PC", seq(nrow(X))))
	return(eig)	
}
