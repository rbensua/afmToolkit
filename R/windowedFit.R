windowedFit <- function(X, width){
  lenX <- dim(X)[1]
  SEQ1 <- seq(width+1,lenX-width)
  SEQ2 <- lapply(SEQ1, function(x) (x-width):(x+width))
  OUT <- lapply(SEQ2, function(a) coef(lm.fit(X[a,1:2],X[a,3]))[2])
  OUT <- base:::simplify2array(OUT,higher=TRUE)
  return(OUT)
}