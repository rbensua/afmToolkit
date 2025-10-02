#' @title afmIndentation
#'
#' @description This function computes the deformation of the sample from the 
#' calibrated Force-Distance curve, by substracting Z to the Zero Force Point 
#' calculated with afmZeroPointSlope function. 
#' @usage afmIndentation(afmdata)
#' @param afmdata An \code{afmdata} object. It should be a valid afmdata
#'  object upon which the Contact Point, the baseline correction and the Zero 
#'  Force Point must have been calculated first (using functions 
#'  \code{afmContactPoint()}, \code{afmBaselineCorrection()}) and
#'   \code{afmZeroPointSlope()}
#'
#' @return Returns a list with one field:
#'
#' \code{afmdata}: An afmdata class in which a \code{Indentation} column is added 
#' in the \code{data} field.
#' @examples 
#' data <- afmReadJPK("force-save-JPK-3h.txt.gz", path = path.package("afmToolkit"))
#' data <- afmContactPoint(data, width = 20, mul1 = 1, mul2 = 20)
#' data <- afmDetachPoint(data, width = 40, mul1 = 3, mul2 = 40)
#' data <- afmBaselineCorrection(data)
#' data <- afmZeroPointSlope(data, segment = "approach")
#' data <- afmIndentation(data)
#' head(data$data)
#' @export


afmIndentation <- function(afmdata){
  if (is.afmexperiment(afmdata)){
    afmdata <- lapply(afmdata, afmIndentation)
    return(afmexperiment(afmdata))
  }else if(is.afmdata(afmdata)){
  if (!("ForceCorrected" %in% names(afmdata$data))) {
    stop("Baseline correction should be done first!")
  }
  if (!any ("Slopes" %in% names(afmdata))) {
    stop("Zero Point Force should be found first (run afmZeroPointSlopes)!")
  }
  if (!afmdata$params$SpringConstant >0) {
    stop("The cantilever spring constant should be provided!")
  }
  Z0Point <- afmdata$Slopes$Z0Point
  CantileverDeflection <- -afmdata$data$ForceCorrected / afmdata$params$SpringConstant
  Indentation <- afmdata$data$Z - Z0Point - CantileverDeflection
  afmdata$data$Indentation <- Indentation
  return(afmdata(afmdata))
  }else{
    stop("Error: input is not a valid afmdata or afmexperiment.")
  }
}