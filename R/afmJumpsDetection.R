jumpsDetection <- function(Z,Force,width=10,Delta=TRUE){
  n = length(Z)
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