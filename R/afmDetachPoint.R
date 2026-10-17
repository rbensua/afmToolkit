#' @title Detach point
#'   
#' @description Find the detach point (or unbinding point) for the Force-Distance curve 
#' following the local regression and two thresholds methods described in Microscopy
#' Research and Technique 2013 (see reference).
#' 
#' The procedure is similar to the one used by the \code{afmContactPoint()} function for
#' obtaining the contact point.
#' @param afmdata A Force-Distance curve with the afmdata structure. It should be a list
#'   with at least the 'data' field with a data frame of at least 4 columns.
#' @param width Width of the window for the local regression (in vector position units)
#' @param mul1 First multiplier for the first alarm threshold
#' @param mul2 Second multiplier for the second alarm threshold
#' @param lagdiff Lag for estimating the differences in Delta (or slopes) signal. By
#'   default it takes the same value as the window with.
#' @param Delta Logical. If TRUE, then the statistic for determining the contact point is
#'   the differences between two consecutive values of the slope of the local regression
#'   line. If FALSE then the slope itself is used.
#' @param loessSmooth Logical If TRUE, a loess smoothing (via loess.smooth()) is done
#'   prior to the determination of the contact point. The span of the  smoothing is 0.05
#'   (5%), the degree is 2 and the number of points equals the number of points in the
#'   approach segment.
#' @param silent Logical. If TRUE supresses a message showing which curve is being processing
#'  (use it when processing a high number of curves). Defaults to FALSE.
#' @return An \code{afmdata} class variable which will consist on the original input
#'   \code{afmdata} variable plus a new list named \code{DP} with the following fields:
#'   
#'   \code{DP} The detach point value.
#'   
#'   \code{iDP} The position in the array for the detach point value.
#'   
#'   \code{delta} The delta signal.
#'   
#'   \code{noise} The noise of the delta signal
#' @examples
#' data <- afmReadJPK("force-save-JPK-3h.txt.gz", path = path.package("afmToolkit"))
#' width <- 10
#' mul1 <- 2
#' mul2 <- 40
#' data <- afmDetachPoint(data, width = width, mul1 = mul1, mul2 = mul2)
#' \dontrun{
#' plot(data, segment = "retract") + geom_vline(xintercept = data$DP$DP, lty = 2)
#' }
#' @references Benitez R., Moreno-Flores S., Bolos V. J. and Toca-Herrera J.L. (2013). "A 
#' new automatic contact point detection algorithm for AFM force curves". Microscopy
#' research and technique, \strong{76} (8), pp. 870-876.
#' @seealso \code{\link{afmContactPoint}}
#' @importFrom stats loess.smooth
#' @export
afmDetachPoint <- function(afmdata,width=1,
                           mul1,
                           mul2, 
                           lagdiff = width, 
                           Delta=TRUE, 
                           loessSmooth = FALSE, 
                           silent = FALSE){
  afmDispatch(
    afmdata,
    afmDetachPointSingle,
    width = width,
    mul1 = mul1,
    mul2 = mul2,
    lagdiff = lagdiff,
    Delta = Delta,
    loessSmooth = loessSmooth,
    silent = silent,
    error = "Error: Input is not a valid afmdata or afmexperiment."
  )
}

# Internal: computes the detach point of a single afmdata curve.
# The afmdata/afmexperiment dispatch is done by afmDetachPoint()
# via the afmDispatch() helper.
afmDetachPointSingle <- function(afmdata,width=1,
                                 mul1,
                                 mul2, 
                                 lagdiff = width, 
                                 Delta=TRUE, 
                                 loessSmooth = FALSE){
  Segment <- NULL
  data.retract <- subset(afmdata$data, Segment == "retract")
  Z <- data.retract$Z
  Force <- data.retract$Force
  n <- length(Z)
  direction <- Z[n] - Z[1]
  if (direction > 0){
    Z <- rev(Z)
    Force <- rev(Force)
  }
  if (loessSmooth){
    data.contact.smoothed <- loess.smooth(Z, Force, span = 0.05, 
                                          degree = 2, evaluation = n)
    Z <- data.contact.smoothed$x
    Force <- data.contact.smoothed$y
    if (Z[n]-Z[1] > 0){
      Z <- rev(Z)
      Force <- rev(Force)
    }
  }
  
  app <- matrix(c(rep(1,n),Z,Force),nrow = n,ncol = 3)
  bRoll <- windowedFit(app,width) 
  delta<- diff(bRoll, lag = lagdiff)
  delta <- c(rep(0,width+lagdiff),delta,rep(0,width))
  if (!Delta){
    delta <- c(rep(0,width),bRoll,rep(0,width))
  }

  noise <- sd(delta[(width):(width+as.integer(0.1*n))], na.rm = TRUE)
  tol1 <- mul1*noise
  tol2 <- mul2*noise
  
  if(tol2 > max(abs(delta))){
    tol2 <- max(abs(delta))-0.05*diff(range(abs(delta)))
  }
  
  idxGrTol2 <- which(abs(delta)>tol2)
  idxSmTol1 <- which(abs(delta)<tol1)
  j <- max(idxSmTol1[idxSmTol1<min(idxGrTol2)])
    
    if ((j > 1) & (delta[j] != 0)){
      eps <- (tol1-abs(delta[j]))/abs(delta[j+1]-delta[j]) # factor de proporcionalidad entre 0 y 1    
    }
    else{
      eps <- 0
    }
    i_detach = min(j+width,n-1)
    z_detach = Z[i_detach]
    z_detach = as.numeric(z_detach+eps*(Z[i_detach+1]-z_detach))
    DP <- list(DP=z_detach, iDP=i_detach,delta=delta,noise=noise)
    return(append.afmdata(afmdata,DP))
}
