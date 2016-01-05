afmOffset <- function(afmdata, offSetLength, units = c("points","meters")){
  n <- nrow(afmdata$data)
  if(units == "points"){
    init <- offSetLength
    end <-  n - offSetLength
  }else {
    firstZ <- afmdata$data$Z[1]
    init <- min(which(abs(afmdata$data$Z-firstZ)>offSetLength))
    end <- max(which(abs(afmdata$data$Z-firstZ)>offSetLength))
  }
  afmdata$data <- afmdata$data[init:end,]
  return(afmdata(afmdata))
}