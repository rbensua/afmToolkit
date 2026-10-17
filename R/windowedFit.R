#' @title Linear fit in a running window
#' 
#' @description 
#' This is an internal function used by the \code{afmContactPoint} 
#' and \code{afmDetachPoint} functions. 
#' It computes the slopes of a linear fit to the data in a window of a given radius.
#' This function should not be used directly unless by experienced users.
#' @usage windowedFit(X, width)
#' @param X Least squares matrix on the form [1 z Force], 
#' according to input parameters in function \code{lm.fit} 
#' @param width Width of the window for the local regression (in vector position units)

#' @return \code{OUT} A vector of length \code{nrow(X)-2*width}, 
#' containing with the slopes of the fits.
#' @examples 
#' n <- 100
#' x <- seq(0,2*pi,length.out = n)
#' y = sin(x)+0.1*rnorm(n)
#' X <- matrix(c(rep(1,n),x,y),nrow = n,ncol = 3)
#' width <- 5
#' b <- windowedFit(X,width)
#' plot(x[(width+1):(n-width)],b,xlab = "x",ylab = "y",type = "l")
#' lines(x,y,col = "red")
#' legend("bottomleft",c("Slopes","Signal"),col = c(1,2),lty = 1)
#' @importFrom stats lm.fit
#' @export
windowedFit <- function(X, width){
  lenX <- dim(X)[1]
  centers <- seq(width + 1, lenX - width)
  nout <- length(centers)
  if (nout == 0) {
    return(numeric(0))
  }
  z <- X[, 2]
  y <- X[, 3]
  m <- 2 * width + 1
  lo <- centers - width
  hi <- centers + width
  # Prefix sums: allows computing every window's fit in O(1),
  # so the whole call is O(n) instead of one lm.fit per window.
  cs_z <- c(0, cumsum(z))
  cs_y <- c(0, cumsum(y))
  cs_zz <- c(0, cumsum(z * z))
  cs_zy <- c(0, cumsum(z * y))
  Sz <- cs_z[hi + 1] - cs_z[lo]
  Sy <- cs_y[hi + 1] - cs_y[lo]
  Szz <- cs_zz[hi + 1] - cs_zz[lo]
  Szy <- cs_zy[hi + 1] - cs_zy[lo]
  # Slope of the least squares fit of y on [1, z] within each window
  OUT <- (m * Szy - Sz * Sy) / (m * Szz - Sz * Sz)
  # lm.fit names the slope coefficient "x2" when the columns of X are
  # unnamed; preserve that behaviour for backward compatibility
  names(OUT) <- rep(if (is.null(colnames(X))) "x2" else colnames(X)[2],
                    nout)
  return(OUT)
}
