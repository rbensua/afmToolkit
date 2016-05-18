#' @title Detach point
#'
#' @description
#' Finds the detach point in for the Force-Distance curve
#' following the method described in Microscopy Research and Technique 2013
#'
#' @param afmdata A Force-Distance curve with the afmdata structure. It should be a list with at least the 'data' field with a data frame of at least 4 columns.
#' @param width Width of the window for the local regression (in vector position units)
#' @param mul1 First multiplier for the first alarm threshold
#' @param mul2 Second multiplier for the second alarm threshold
#' @param Delta Logical. If TRUE, then the statistic for determining the contact point is the
#' differences between two consecutive values of the slope of the local regression line.
#' If FALSE then the slope itself is used.
#' @param loessSmooth Logical If TRUE (default), a loess smoothing (via loess.smooth()) is done prior to
#' the determination of the contact point. The span of the  smoothing is 0.05 (5%), the degree is 2 and the
#' number of points equals the number of points in the approach segment.
#' @return A list of:
#'
#' \code{CP} The contact point value.
#'
#' \code{iCP} The position in the array for the contact point value.
#'
#' \code{delta} The delta signal.
#'
#' \code{noise} The noise of the delta signal
#' @examples
#' JPKexample <- paste(path.package("afmToolkit"), "force-save-JPK-3h.txt",sep = "/")
#' data <- afmReadJPK(JPKexample)
#' width <- 20
#' mul1 <- 1
#' mul2 <- 10
#' CP <- afmContactPoint(data, width = width, mul1 = mul1, mul2 = mul2)
#' plot(data, segment = "approach") + geom_vline(xintercept = CP$CP, lty = 2)
#' @export
afmDetachPoint <- function(afmdata,width=1,mul1,mul2, lagdiff = width, 
                            Delta=TRUE, loessSmooth = TRUE){
  if (is.afmexperiment(afmdata)){
    DP <- lapply(afmdata, function(x) afmDetachPoint(x, width = width,
                                                      mul1 = mul1,
                                                      mul2 = mul2,
                                                      lagdiff = lagdiff,
                                                      Delta = Delta,
                                                      loessSmooth = loessSmooth))
    afmexperiment <- mapply(afmdata,DP, FUN = function(x,y) append.afmdata(x,y, name = "DP"),
                            SIMPLIFY = FALSE)
    return(afmexperiment(afmexperiment))
  }else if (is.afmdata(afmdata)){
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
  
  
  b<- array(0,dim=c(n,1))
  delta <- array(0,dim=c(n,1))
  imax <- n-width
  app <- matrix(c(rep(1,n),Z,Force),nrow = n,ncol = 3)
  # bRoll <- rollapply(app,width,FUN = function(X) 
  #  coef(lm.fit(X[,1:2],X[,3]))[2],by.column=FALSE,align="right")
  bRoll <- windowedFit(app,width) 
  delta<- diff(bRoll, lag = lagdiff)
  delta <- c(rep(0,width+lagdiff),delta,rep(0,width))
  if (!Delta){
    delta <- c(rep(0,width),bRoll,rep(0,width))
  }
  
#  noise <- sd(delta[(width+as.integer(n/3)):(width+as.integer(n/3)+as.integer(0.1*n))],na.rm=TRUE)
  noise <- sd(delta[(width):(width+as.integer(0.1*n))], na.rm = TRUE)
  tol1 <- mul1*noise
  tol2 <- mul2*noise
  
  if(tol2 > max(abs(delta))){
    tol2 <- max(abs(delta))-0.05*diff(range(abs(delta)))
  }
  
  idxGrTol2 <- which(abs(delta)>tol2)
  idxSmTol1 <- which(abs(delta)<tol1)
  if (length(idxGrTol2>tol2) ==0){
    cat("mul2 is too large. Set a smaller value for mul2\n")
    return(list(CP=NA,iCP=NA,delta=delta,noise=noise))
    
  } else{
    j <- max(idxSmTol1[idxSmTol1<min(idxGrTol2)])#+1  
    
    if ((j > 1) & (delta[j] != 0)){
      eps <- (tol1-abs(delta[j]))/abs(delta[j+1]-delta[j]) # factor de proporcionalidad entre 0 y 1    
    }
    else{
      eps <- 0
    }
    i_detach = min(j+width,n-1)
    z_detach = Z[i_detach]
    z_detach = as.numeric(z_detach+eps*(Z[i_detach+1]-z_detach));  
    return(list(DP=z_detach, iDP=i_detach,delta=delta,noise=noise))
  }
  }else{
    stop("Error: Input is not a valid afmdata or afmexperiment.")
  }
}