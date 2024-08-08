#' Color opacity.
#' 
#' @description
#' Makes colors opaque / transparent.
#' 
#' @param color a vector of colors (color names, hexadecimal strings or integers from [grDevices::palette]).
#' @param alpha a numeric vector, degree of opacity (from 0 to 1), recycled if necessary.
#' @param returns A vector of modified colors.
#' @param details This function defines opacity of a color, 0 being fully transparent and 1 being fully visable.
#' @export

opacity <- function(color, alpha) {
	hex <- c("0","1","2","3","4","5","6","7","8","9","A","B","C","D","E","F")
	num2hex <- function(x, hex) paste0(hex[(x - x %% 16) / 16 + 1], hex[x %% 16 + 1])
	if (length(color) != length(alpha) & !any(c(length(color), length(alpha)) == 1)) {
		stop("Vector lengths not correct")
	} else if (length(color) == 1 & length(alpha) > 1) {
		color <- rep(color, length(alpha))
	} else if (length(alpha) == 1 & length(color) > 1) {
		alpha <- rep(alpha, length(color))
	}
	alpha <- as.integer(255 * alpha)
	rgb <- rbind(grDevices::col2rgb(color), alpha)
	alphacolor <- paste0("#", apply(apply(rgb, 2, num2hex, hex=hex), 2, paste, collapse=""))
	return(alphacolor)
}



