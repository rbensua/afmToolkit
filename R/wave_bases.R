wave_bases <- function(mother = c('MORLET','PAUL','DOG'),k,scale,param=NULL){
  mother <- toupper(mother)
  mother <- match.arg(mother)
  n <- length(k)
  
  if(mother == "MORLET"){
    if (is.null(param)){ 
      param <- 6
    }
    k0 <- param
    expnt <- -(scale*k-k0)^2/2*(k>0)
    Norm <- sqrt(scale*k[2])*pi^(-.25)*sqrt(n)
    daughter <- Norm*exp(expnt)
    daughter <- daughter*(k>0)
    fourier_factor <- (4*pi)/(k0+sqrt(2+k0^2))
    coi <- fourier_factor/sqrt(2)
  }else if(mother == "PAUL"){
    if (is.null(param)){ 
      param <- 4
    }
    m <- param
    expnt <- -(scale*k)*(k>0)
    Norm <- sqrt(scale*k[2])*(2^m/sqrt(m*prod(2:(2*m-1))))*sqrt(n)
    daughter <- Norm*((scale*k)^m)*exp(expnt)
    daughter <- daughter*(k>0)
    fourier_factor <- 4*pi/(2*m+1)
    coi  <- fourier_factor*sqrt(2)
  }else if(mother == "DOG"){
    if (is.null(param)){ 
      param <- 2
    }
    m <- param
    expnt <- -(scale*k)^2/2
    Norm <- sqrt(scale*k[2]/gamma(m+0.5))*sqrt(n)
    daughter <- -Norm*(1i^m)*((scale*k)^m)*exp(expnt)
    fourier_factor <- 2*pi*sqrt(2/(2*m+1))
    coi <- fourier_factor/sqrt(2)
  }
  return(list(daughter = daughter, fourier_factor = fourier_factor, coi = coi))
}