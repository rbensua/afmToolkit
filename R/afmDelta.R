afmDelta <- function(afmdata,  width = 1, segment =
                       c("approach","retract","both"), 
                     lagdiff = width) {
  segment <- match.arg(segment)
  data <- subset(afmdata$data, Segment == segment)
  Z <- data$Z
  if (!("ForceCorrected" %in% names(afmdata$data))) {
    Force <- data$Force
  } else{
    Force <- data$ForceCorrected
  }
  n = length(Z)
  delta <- array(0,dim = c(n,1))
  app <- matrix(c(rep(1,n),Z,Force),nrow = n,ncol = 3)
  bRoll <- windowedFit(app,width)
  delta <- diff(bRoll, lag = lagdiff)
  delta <- c(rep(0,width + lagdiff),delta,rep(0,width))
  bRoll <- c(rep(0,width),bRoll,rep(0,width))
  
  return(list(delta = delta, bRoll = bRoll))
}