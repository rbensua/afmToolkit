core <- function(x,n){
  if (n>length(x)){
    stop("Error. The length of the vector x should be greater than n")
  }
  if (n == length(x)){
   extracted <- x 
  }else{
  n1 <- floor((length(x)-n)/2)
  n2 <- length(x)-n1
  extracted <- x[(n1+1):(n2-1)]
  }
  return(extracted)
}