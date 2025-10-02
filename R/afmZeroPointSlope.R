#' @title Zero Force Point and Slope
#'   
#' @description This function finds the point of zero force (real contact point) and the
#'   slope of the contact part of the Force-Distance curve.
#' @param afmdata An \code{afmdata} object. It should be a valid afmdata object upon which
#'   the Contact Point and the baseline correction must have been calculated first (using
#'   functions \code{afmContactPoint()} and \code{afmBaselineCorrection()})
#' @param fstar Value such that fstar * sd is to be considered as zero Force, where sd is
#'   the standard deviation of Force at the basline. It takes fstar = 0 as default value,
#'   meaning that zero force is actually zero.
#' @param segment The segment on which everything is calculated.
#' @param silent Logical value. If TRUE suppresses the message indicating the name of the curve being processed
#' (useful for batch-processing large number of curves). Defaults to FALSE
#'   
#' @return An \code{afmdata} class variable which will consist on the original input
#'   \code{afmdata} variable plus a new list named \code{Slopes} with the following
#'   fields: \code{Z0Point}: Point of zero force. \code{Slope}: Slope of the best fit line
#'   in the contact part of the Force-Distance curve.
#' @examples
#' data <- afmReadJPK("force-save-JPK-2h.txt.gz", path = path.package("afmToolkit"))
#' data <- afmContactPoint(data, width = 20, mul1 = 1, mul2 = 20)
#' data <- afmDetachPoint(data, width = 40, mul1 = 3, mul2 = 40)
#' data <- afmBaselineCorrection(data)
#' data <- afmZeroPointSlope(data, segment = "approach")
#' \dontrun{
#' plot(data, segment = "approach") + geom_vline(xintercept = data$Slopes$Z0Point, lty = 2)
#' }
#' @export
afmZeroPointSlope <-
  function(afmdata, fstar = 0, segment = c("approach", "retract"), silent = FALSE) {
    Segment <- NULL
    if (is.afmexperiment(afmdata)) {
      afmdata <-
        lapply(afmdata, function(x){
          if(!is.null(x$params$curvename)& !silent){
            print(paste("Processing curve: ",x$params$curvename), sep = " ")
          }
          afmZeroPointSlope(x, fstar = fstar, segment = segment)
        })
      return(afmexperiment(afmdata))
    } else if (is.afmdata(afmdata)) {
      if (!("ForceCorrected" %in% names(afmdata$data))) {
        stop("Baseline correction should be done first!")
      }
      # if (!any(sapply(afmdata, function(x)
      #   "CP" %in% names(x)))) {
      if(!"CP" %in% names(afmdata)){
        stop("Contact Point should be found first!")
      }
      segment <- match.arg(segment)
    if(segment == "approach" & is.afmmulti(afmdata)){
      segment <- "approach1"
    }
      Z <- subset(afmdata$data, Segment == segment)$Z
      isTime <- "Time" %in% colnames(afmdata$data)
      
      if (isTime) {
        Time <- subset(afmdata$data, Segment == segment)$Time  
      }
      ForceCorrected <-
        subset(afmdata$data, Segment == segment)$ForceCorrected
      Zmin <- Z[which.min(ForceCorrected)]
      Fmin <- min(ForceCorrected)
      if (grepl("approach",segment)) {
        # Computing the standard deviation in the baseline
        stddev <- sd(ForceCorrected[which(Z > afmdata$CP$CP)])
        zerovalue <- fstar*stddev;
        #indicesSlope <-
        #  which(ForceCorrected > 0 & Z < min(afmdata$CP$CP, Zmin))
        indicesSlope <-
          which(ForceCorrected <= zerovalue & Z <=  Zmin)
        
        i1 <- max(indicesSlope)
        #i0 <- i1 - 1
        i0 <- i1+1
        # if (abs(Zmin) < abs(afmdata$CP$CP)) {
        #   Z0Point <- Z[i0] - ForceCorrected[i0] * (Z[i1] - Z[i0]) /
        #     (ForceCorrected[i1] - ForceCorrected[i0])
        # } else {
        #   Z0Point <- afmdata$CP$CP
        # }
        Z0Point <- Z[i0] + (zerovalue - ForceCorrected[i0]) * (Z[i1] - Z[i0]) /
          (ForceCorrected[i1] - ForceCorrected[i0])
        if (isTime){
          t0Point <- Time[i0] + (zerovalue - ForceCorrected[i0]) * (Time[i1] - Time[i0]) /
            (ForceCorrected[i1] - ForceCorrected[i0])
        } else{
          t0Point <- NULL
        }
        
      } else {
        indicesSlope <-
          which(ForceCorrected > 0 & Z < min(afmdata$DP$DP, Zmin))
        i1 <- max(indicesSlope)
        i0 <- i1 + 1
        #  if (Zmin < afmdata$DP$DP){
        Z0Point <- Z[i0] - ForceCorrected[i0] * (Z[i1] - Z[i0]) /
          (ForceCorrected[i1] - ForceCorrected[i0])
        if (isTime){
          t0Point <- Time[i0] - ForceCorrected[i0] * (Time[i1] - Time[i0]) /
            (ForceCorrected[i1] - ForceCorrected[i0])
        } else{
          t0Point <- NULL
        }
        
        # } else {
        #   Z0Point <- afmdata$DP$DP
        #  }
      }
      indicesSlope2 <-
        which(Z <=  Zmin)
      Zslope <- Z[indicesSlope2]
      Fslope <- ForceCorrected[indicesSlope2]
      FitSlope <- lm(Fslope ~ Zslope)
      slope <- coef(FitSlope)[2] # Second coefficient of the fit
      Slopes <- list(Z0Point = Z0Point, t0Point = t0Point, Slope = slope)
      return(append.afmdata(afmdata,Slopes))
    } else{
      stop("Error: input is not a valid afmdata or afmexperiment.")
    }
  }
