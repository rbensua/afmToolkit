#' @title AFM experiment
#'
#' @description
#' This function creates an \code{afmexperiment} structure, which is as list 
#' (or an array) of elements of \code{afmdata} class.
#'
#' @param \code{data}: A variable of  \code{afmdata} class, or a list of 
#' elements of \code{afmdata} class.
#' @param \code{ID}: Character string with the identifier of 
#' the \code{data} variable or a string array in case \code{data} is a 
#' list of \code{afmdata} variables. 
#' 
#' @return An object of class \code{afmexp}.
#'
#' @examples
#'#Making some artifical data following a L-J 12-6 potential
#'n <- 1000
#'z <- seq(from = 9e-3, to = 1e-1, length.out = n )
#'u0 <- 1e-5
#'z0 <- 1e-2
#'Force <- -u0*(12*z0^6/z^7-12*z0^12/z^13)
#'Segment <- rep("approach",n)
#'AFMcurve <- afmdata(data.frame(Z = z, Force = Force, Segment  = Segment))
#'plot(AFMcurve)


afmexperiment <-
  function(data, ID=NULL) {
    if (!is.afmexperiment(data)) {
        if (is.afmdata(data)){
          experiment <- list(data)
          if (is.null(ID)){
           names(experiment) <- "1" 
          }else {
          names(experiment) <- ID
          }
          return(structure(experiment, class = "afmexperiment"))
        }else if (is.list(data)){
          if (prod(sapply(data, is.afmdata))){
            experiment <- data
            if(is.null(ID)){
              if (is.null(names(data))){
                names(experiment) <- as.character(seq_along(data))
              }
            }else{
              names(experiment) <- ID
            }
          }else{
            stop("All elements in the list should be of afmdata class.")
          }
          return(structure(experiment, class = "afmexperiment"))
        }      
    } else{
      return(data)
    }
    
}