#'@title AFM data
#'  
#'@description This function creates an \code{afmdata} structure, which is as list with at
#'least one field called \code{data} which is a data frame with a valid AFM data, that is,
#'at least 3 variables called "Z", "Force", and "Segment".
#'
#'@usage afmdata(data, dstr = "Z", Fstr = "Force", Segstr = "Segment", tstr = "Time", 
#'  params = list(SpringConstant = numeric(), curvename = NULL ))
#'@param data A data frame consisting in 3 or 4 columns. A minimum of "Z" (or "distance"),
#'  "Force" and "Segment". Optionally a fourth column with "Time" could be added.
#'@param dstr Character string with the posible names for the distance variable.
#'@param Fstr Character string with the posible names for the force variable.
#'@param Segstr Character string with the posible names for the Segment variable.
#'@param tstr Character string with the posible names for the time variable.
#'@param params A list that may contain parameters describing the F-d curve. At least will
#'  contain the \code{SpringConstant} and the \code{curvename}, being the former the
#'  cantilever spring constant and the latter  a F-d curve ID. Function \code{afmReadJPK}
#'  will try to obtain the spring constant from the file header and the curvename from the
#'  data file name.
#'@return An object of class \code{afmdata}
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
#'@seealso \code{\link{afmexperiment}}
#'@export

afmdata <-
  function(data, dstr = "Z", Fstr = "Force", Segstr = "Segment", tstr = "Time",
           params = list(SpringConstant = numeric(), curvename = NULL )) {
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
