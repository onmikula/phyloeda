pcoord <- function(D) {
	n <- nrow(D)
	C <- diag(n) - (1 / n) * rep(1, n) %*% t(rep(1, n))
	A <- -0.5 * D^2
	delta <- C %*% A %*% C
	eig <- eigen(delta)
	pos <- eig$values >= 0
	val <- eig$values[pos]
	vec <- eig$vectors[,pos] %*% diag(sqrt(val))
	dimnames(vec) <- list(rownames(D), paste0("PCo", seq(sum(pos))))
	result <- list(vectors=vec, values=val)
	class(result) <- "pcoa"
	return(result)
}
