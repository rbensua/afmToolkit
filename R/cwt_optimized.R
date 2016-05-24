cwt_optimized <- function(scales, signal, wname = c("MORLET","HAAR")){
  wname <- toupper(wname)
  precision <- max(12,floor(log2(max(scales)))+4)
  n <- length(signal)
  ns <- length(scales)
  coefs <- matrix(0,nrow = ns, ncol = n)
  if (wname == "MORLET"){
    coefs <- coefs + 1i*coefs
    k <- 1:trunc(n/2)
    k <- k*2*pi/n
    k <- c(0,k,-k[trunc((n-1)/2):1])
    f <- fft(signal)
    for (i in 1:ns){
      wav <- wave_bases("Morlet",k,scales[i])
      coefs[i,] <- fft(f*wav$daughter, inverse = TRUE)/length(f)
    }
  } else if (wname == "HAAR"){
    step <- 2^(-precision)
    xval <- seq(from = 0, to =1+step, by = step)
    xmax <- max(xval)
    psi_integ <- xval*(xval<1/2)+(1-xval)*(xval>=1/2)
    for (k in 1:ns){
      a <- scales[k]
      j <- 1+floor((0:a*xmax)/(a*step))
      if (length(j) == 1){
        j <- c(1,1)
      } 
      f <- rev(psi_integ[j])
      coefs[k,] <- -sqrt(a)*core(diff(convolve(signal,f,type = "open")),n)
    }
  }
  return(coefs)
}