#' @title Zero Force Point and Slope
#'
#' @description This function finds the point of zero force (real contact point)
#' and the slope of the contact part of the Force-Distance curve.
#'
#' @param \code{afmdata}: An \code{afmdata} object. It should be a valid afmdata object upon
#' which the Contact Point and the baseline correction must have been calculated first
#' (using functions \code{afmContactPoint()} and \code{afmBaselineCorrection()})
#' @param \code{segment}: The segment on which everything is calculated.
#'
#' @return Returns a list with two fields:
#'
#' \code{Z0Point}: Point of zero force.
#' \code{Slope}: Slope of the best fit line in the contact part of the Force-Distance curve.
#' @export





afmZeroPointSlope <- function(afmdata, segment = c("approach","retract")) {
  if (!("ForceCorrected" %in% names(afmdata$data))) {
    stop("Baseline correction should be done first!")
  }
  if (!any(sapply(afmdata, function(x)
    "CP" %in% names(x)))) {
    stop("Contact Point should be found first!")
  }
  segment <- match.arg(segment)
  Z <- subset(afmdata$data, Segment == segment)$Z
  ForceCorrected <-
    subset(afmdata$data, Segment == segment)$ForceCorrected
  Zmin <- Z[which.min(ForceCorrected)]
  Fmin <- min(ForceCorrected)
  
  if (segment == "approach") {
    indicesSlope <-
      which(ForceCorrected > 0 & Z < min(afmdata$CP$CP, Zmin))
    i1 <- min(indicesSlope)
    i0 <- i1 - 1
    if (abs(Zmin) < abs(afmdata$CP$CP)){
      Z0Point <- Z[i0] - ForceCorrected[i0] * (Z[i1] - Z[i0]) /
        (ForceCorrected[i1] - ForceCorrected[i0])
    } else {
      Z0Point <- afmdata$CP$CP
    }
    
  } else {
    indicesSlope <-
      which(ForceCorrected > 0 & Z < min(afmdata$DP$DP, Zmin))
    i1 <- max(indicesSlope)
    i0 <- i1 + 1
  #  if (Zmin < afmdata$DP$DP){
      Z0Point <- Z[i0] - ForceCorrected[i0] * (Z[i1] - Z[i0]) /
        (ForceCorrected[i1] - ForceCorrected[i0])
   # } else {
   #   Z0Point <- afmdata$DP$DP
  #  }
  }
  
  Zslope <- Z[indicesSlope]
  Fslope <- ForceCorrected[indicesSlope]
  FitSlope <- lm(Fslope ~ Zslope)
  slope <- coef(FitSlope)[2] # Second coefficient of the fit
  return(list = list(Z0Point = Z0Point, Slope = slope))
}