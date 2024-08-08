#' Elbow method.
#' 
#' @description
#' An elbow method to determine the number of important principal components.
#' 
#' @param y numeric, eigenvalues from PCA.
#' @param x ranks of the principal components (defaults to `seq_along(y)`).
#' @returns A list with components `x` and `y` (the input data), `b` (breakpoints), `rmse` (root mean squared error)
#'   and `best` (the optimal breakpoint, corresponding to `b[which.min(rmse)]`). It is assigned a class `elbow`.
#' @details The method finds a breakpoint separating first principal components accounting for most of the variation
#'   from the remaining ones. The function is an implementation of approach described by Salvador and Chan 2004
#'   (https://doi.org/10.1109/ICTAI.2004.50). It can be used to find breakpoints also in other analogous data.
#' @export

elbow <- function(y, x=seq_along(y)) {
	k <- length(y)
	b <- 2:(k - 2)
	w <- rbind(b, (k-b)) / k
	lms <- lapply(b, function(p) list(L=lm(y[1:p] ~ x[1:p]), R=lm(y[(p+1):k] ~ x[(p+1):k])))
	rmse <- colSums(w * sqrt(sapply(lms, function(m) sapply(m, deviance))))
	result <- list(x=x, y=y, b=b, rmse=rmse, best=b[which.min(rmse)])
	class(result) <- "elbow"
	return(result)
}


#' Elbow method.
#' 
#' @description
#' A plot method for the result of [elbow].
#' 
#' @param x an object of class `elbow`, result of [elbow].
#' @export

plot.elbow <- function(x) {
	par(mai=c(1.02, 1.02, 0.42, 0.42))
	plot(x$x, x$y, ylim=range(c(x$y, x$rmse)), xlab="breakpoints", ylab="RMSE", type="b", col="darkgrey", cex.lab=1.5, cex.axis=1.15, pch=21, cex=1.15, lwd=1.5)
	lines(x$x[2:(length(x$x)-2)], x$rmse, type="b", pch=21, col="black", bg=ifelse(seq_along(x$y) == (x$best-1), "black", "white"), cex=1.15, lwd=1.5)
}
