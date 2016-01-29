
afmDetachPoint <- function(afmdata,width=1,mul1,mul2, lagdiff = width, 
                            Delta=TRUE, loessSmooth = TRUE){
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
}