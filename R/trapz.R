trapz <- function(x,y){
  idx = 2:length(x)
  return(as.double((x[idx] - x[idx - 1]) %*% 
                     (y[idx] + y[idx - 1]))/2)
}