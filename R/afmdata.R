#' @title AFM data
#'
#' @description
#' This function creates an \code{afmdata} structure, which is as list with at
#'  least one field called \code{data} which is a data frame with a valid AFM 
#'  data, that is, at least 3 variables called "Z", "Force", and "Segment".
#'
#' @param \code{data}: A data frame consisting in 3 or 4 columns. A minimum of 
#' "Z" (or "distance"), "Force" and "Segment". Optionally a fourth column with 
#' "Time" could be added.
#' @param \code{dstr}: Character string with the posible names for the distance
#'  variable.
#' @param \code{Fstr}: Character string with the posible names for the force 
#' variable.
#' @param \code{Segstr}: Character string with the posible names for the Segment
#'  variable.
#' @param \code{tstr}: Character string with the posible names for the time
#'  variable.
#'
#' @return An object of class \code{afmdata}
#'
#'@examples
#'#Making some artifical data following a L-J 12-6 potential
#'n <- 1000
#'z <- seq(from = 9e-3, to = 1e-1, length.out = n )
#'u0 <- 1e-5
#'z0 <- 1e-2
#'Force <- -u0*(12*z0^6/z^7-12*z0^12/z^13)
#'Segment <- rep("approach",n)
#'AFMcurve <- afmdata(data.frame(Z = z, Force = Force, Segment  = Segment))
#'plot(AFMcurve)

afmdata <-
  function(data, dstr = "Z", Fstr = "Force", Segstr = "Segment", 
           params = list(SpringConstant = numeric())) {
    if (!is.afmdata(data)) {
      if (is.data.frame(data)) {
        if (!prod(c(dstr,Fstr,Segstr) %in% names(data))) {
          stop("data should be a valid AFM data frame")
        } else {
          return(structure(list(data = data, params = params),class = "afmdata"))
        }
      } else {
        if (!is.list(data)) {
          stop("Not valid data (not a list)")
        } else {
          if (!"data" %in% names(data)) {
            stop("Input does not contain a 'data' field")
          } else {
            if (!is.data.frame(data$data)) {
              stop("Input 'data' field is not a data.frame")
            } else {
              if (!prod(c(dstr,Fstr,Segstr) %in% names(data$data))) {
                stop("'data' field is not a valid AFM data frame")
              } else {
                if (!"params" %in% names(data)){
                return(structure(append(data, list(params = params)),class = "afmdata"))
                } else{
                  return(structure(data,class = "afmdata"))
                }
              }
            }
          }
        }
      }
      
    } else {
      return(data)
    }
  }
