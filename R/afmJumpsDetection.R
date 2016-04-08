
afmJumpsDetection <- function(afmdata,width=10,Delta=TRUE, lessSmooth = TRUE){
  
  data.retract <- subset(afmdata$data, Segment == "retract")
  n <- nrow(data.retract)
  direction <- data.retract$Z[n] - data.retract$Z[1]
  if (loessSmooth){
    data.retract.smoothed <- loess.smooth(data.retract$Z, data.retract$Force,
                                           span = 0.05, degree = 2, evaluation = n)
    Z <- data.retract.smoothed$x
    Force <- data.retract.smoothed$y
    if (direction <0){
      Z <- rev(Z)
      Force <- rev(Force)
    }
  } else{
    Z <- retract$Z
    Force <- retract$ForceCorrected
  }
  if(length(Force)!=n){
    stop("Z and Force must be the same length")
  }
  
  b<- array(0,dim=c(n,1))
  delta <- array(0,dim=c(n,1))
  imax <- n-width
  app <- matrix(c(rep(1,n),Z,Force),nrow = n,ncol = 3)

  bRoll <- windowedFit(app,width) 
  delta<- diff(bRoll)
  delta <- c(rep(0,width+1),delta,rep(0,width))
  if (!Delta){
    delta <- c(rep(0,width),bRoll,rep(0,width))
  }
  return(delta)
}