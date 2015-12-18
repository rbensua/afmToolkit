<<<<<<< HEAD
#' @title Contact point
#' 
#' @description 
=======
>>>>>>> 05e3141ced9de0ca59a9494352758e57a829facc
#' Find the contact point in for the Force-Distance curve
#' following the method described in Microscopy Research and Technique 2013
#' 
#' @param afmdata A Force-Distance curve with the afmdata structure. It should be a list with at least the 'data' field with a data frame of at least 4 columns.
#' @param width Width of the window for the local regression (in vector position units)
#' @param mul1 First multiplier for the first alarm threshold
#' @param mul2 Second multiplier for the second alarm threshold  
#' @param Delta Logical. If TRUE, then the statistic for determining the contact point is the differences between two consecutive values of the slope of the local regression line. If FALSE then the slope itself is used.
<<<<<<< HEAD
#' @return A list of: 
#' 
#' \code{CP} The contact point value. 
#' 
#' \code{iCP} The position in the array for the contact point value. 
#' 
#' \code{delta} The delta signal. 
#' 
#' \code{noise} The noise of the delta signal
#' @export



afmContactPoint <- function(afmdata,width=1,mul1,mul2, lagdiff = width, 
                            Delta=TRUE){
  data.approach <- subset(afmdata$data, Segment == "approach")
  Z <- data.approach$Z
  Force <- data.approach$Force
=======
#' @return A list of: CP The contact point value. iCP The position in the array for the contact point value. delta The delta signal. noise The noise of the delta signal



afmContactPoint <- function(afmdata,width=1,mul1,mul2,Delta=TRUE){
  data.approach <- subset(afmdata$data, Segment == "approach")
  Z <- data.approach$Z
  Force <- data.approach$F
>>>>>>> 05e3141ced9de0ca59a9494352758e57a829facc
  n = length(Z)
  if(length(Force)!=n){
    stop("Z and Force must be the same length")
  }
  b<- array(0,dim=c(n,1))
  delta <- array(0,dim=c(n,1))
  imax <- n-width
  app <- matrix(c(rep(1,n),Z,Force),nrow = n,ncol = 3)
 # bRoll <- rollapply(app,width,FUN = function(X) 
  #  coef(lm.fit(X[,1:2],X[,3]))[2],by.column=FALSE,align="right")
 bRoll <- windowedFit(app,width) 
<<<<<<< HEAD
 delta<- diff(bRoll, lag = lagdiff)
 delta <- c(rep(0,width+lagdiff),delta,rep(0,width))
=======
 delta<- diff(bRoll)
 delta <- c(rep(0,width+1),delta,rep(0,width))
>>>>>>> 05e3141ced9de0ca59a9494352758e57a829facc
  if (!Delta){
    delta <- c(rep(0,width),bRoll,rep(0,width))
  }
 
  noise <- sd(delta[(width+as.integer(n/3)):(width+as.integer(n/3)+as.integer(0.1*n))],na.rm=TRUE)
  tol1 <- mul1*noise
  tol2 <- mul2*noise
  
  idxGrTol2 <- which(abs(delta)>tol2)
  idxSmTol1 <- which(abs(delta)<tol1)
  if (length(idxGrTol2>tol2) ==0){
    cat("mul2 is too large. Set a smaller value for mul2\n")
    return(list(CP=NA,iCP=NA,delta=delta,noise=noise))

  } else{
<<<<<<< HEAD
    j <- max(idxSmTol1[idxSmTol1<min(idxGrTol2)])#+1  
=======
    j <- max(idxSmTol1[idxSmTol1<min(idxGrTol2)])+1  
>>>>>>> 05e3141ced9de0ca59a9494352758e57a829facc
  
  if ((j > 1) & (delta[j] != 0)){
    eps <- (tol1-abs(delta[j]))/abs(delta[j+1]-delta[j]) # factor de proporcionalidad entre 0 y 1    
  }
  else{
    eps <- 0
  }
  i_contact = min(j+width,n-1)
  z_contact = Z[i_contact]
  z_contact = z_contact+eps*(Z[i_contact+1]-z_contact);  
  return(list(CP=z_contact,iCP=i_contact,delta=delta,noise=noise))
}
}