#' @title afmYoungModulus
#'
#' @description This function computes the Young's Modulus of the sample 
#' from the approach curve using different contact models and for different 
#' tip geometries. 
#'
#' @param \code{afmdata}: An \code{afmdata} object. It should be a valid afmdata
#'  object upon which the Contact Point, the baseline correction and the Zero 
#'  Force Point must have been calculated first (using functions 
#'  \code{afmContactPoint()}, \code{afmBaselineCorrection()}) and
#'   \code{afmZeroPointSlope()}
#'
#' @return Returns a list with one field:
#'
#' \code{afmdata}: An afmdata class in which a Deformation vector is added 
#' in the \code{data} field
#' @export


afmYoungModulus <- function(afmdata, thickness = NULL, model = "Hertz", geometry = "pyramid",silent = TRUE, params){
  if (!("ForceCorrected" %in% names(afmdata$data))) {
    stop("Baseline correction should be done first!")
  }
  if (!any ("Slopes" %in% names(afmdata))) {
    stop("Zero Point Force should be found first (run afmSlopes)!")
  }
  if (!afmdata$params$SpringConstant >0) {
    stop("The cantilever spring constant should be provided!")
  }
  if (!("Indentation" %in% names(afmdata$data))){
    stop("Indentation should be computed first(run afmIndentation)!")
  }
  if (is.null(thickness)){
  fitdata <- subset(afmdata$data, Segment == "approach" & Indentation <0)
  } else{
    fitdata <- subset(afmdata$data, Segment == "approach" & Indentation <0 & 
                        abs(Indentation)< thickness)
  }
  fitYM <- lm(ForceCorrected~I(Indentation^2)-1, data = fitdata)
  slope <- coef(fitYM)
  if (!silent){
    print(summary(fitYM))}
  YoungModulus <- as.numeric(slope * sqrt(2)*(1-params$mu^2)/tan(params$alpha*pi/180))
  afmdata$params$YoungModulus <- YoungModulus
  return(afmdata(afmdata))
}