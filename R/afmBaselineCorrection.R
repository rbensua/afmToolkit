#' @title Performs a baseline correction to an AFM F-z curve
#'
#' @description
#' This function performs the baseline correction to an AFM F-z curve within an
#' \code{afmdata} structure.
#'
#'
#' It substracts a best fit line to the curve: for the approach segment,
#' it fits a line to the approach curve points for which |z|>ZPointApp,
#' and for the retract segment it fits a line to the retract curve where
#' |z|>ZPointRet. Every other segment of the curve (e.g. "pause" or
#' "modulation") is corrected according to the \code{fitpause} option:
#' corrected with the approach fit, with the retract fit, or left
#' uncorrected ("none"). When the curve has no retract segment, only
#' the approach segment is fitted and every other segment follows
#' \code{fitpause}.
#'
#' @param afmdata An \code{afmdata} structure.
#' @param ZPointApp Point in the approach segment of the curve
#' that defines the approach baseline
#' @param ZPointRet Point in the retract segment of the curves that 
#' defines the retract baseline
#' #' @param fitpause Behaviour for the baseline correction of every segment
#' other than the approach and retract segments (e.g. "pause",
#' "modulation"): if "approach" (default), those segments are corrected
#' using the best line fit done on the approach segment, if "retract" the
#' best line fit of the retract segment is used, if "none", no baseline
#' correction is done on them.
#' @param vsTime Logical. If TRUE then the baseline correction is performed following the Force vs time approach 
#' described by S. Moreno-Flores (\cite{Moreno Flores (2016)}). 
#' @param sinusoidal Logical. If TRUE, then a pure harmonic is added to the fit of the baseline in order to remove possible periodical noise. Defaults to FALSE.
#' @return \code{afmdata} An \code{afmdata} structure identical to the one in 
#' the input, but with an additional \code{ForceCorrected} column in the
#'  \code{data} dataframe of the \code{afmdata} structure.
#' @importFrom stats lm predict nls
#' @examples
#' AFMcurve <- afmReadJPK("force-save-JPK-2h.txt.gz", path = path.package("afmToolkit"))
#' ZPointApp <- 6.43e-6
#' ZPointRet <- 6.45e-6
#' AFMcurve <- afmBaselineCorrection(AFMcurve,ZPointApp = ZPointApp,ZPointRet = ZPointRet)
#' plot(AFMcurve)
#'
#' # Without providing ZPointApp
#' AFMcurve <- afmReadJPK("force-save-JPK-3h.txt.gz", path = path.package("afmToolkit"))
#' AFMcurve <- afmContactPoint(AFMcurve,width = 10,mul1 = 1,mul2 = 20, 
#'                              loessSmooth = FALSE)
#' AFMcurve <- afmBaselineCorrection(AFMcurve)
#' plot(AFMcurve)
#' 
#' @section References: 
#' Moreno Flores (2016). 
#' Baseline correction of AFM force curves in the force-time representation.
#' Microscopy Research and Technique, 79, (11), pp. 1045-1049.
#' 
#' @export
afmBaselineCorrection <-
  function(afmdata,
           ZPointApp = NULL,
           ZPointRet = NULL, 
           fitpause = c("approach","retract","none"), 
           vsTime = FALSE,
           sinusoidal = FALSE) {
    fitpause <- match.arg(fitpause)
    afmDispatch(
      afmdata,
      afmBaselineCorrectionSingle,
      ZPointApp = ZPointApp,
      ZPointRet = ZPointRet,
      fitpause = fitpause,
      vsTime = vsTime,
      sinusoidal = sinusoidal,
      silent = TRUE,
      error = "input is not a valid afmdata or afmexperiment."
    )
  }

# Internal: baseline correction of a single afmdata curve.
# The afmdata/afmexperiment dispatch is done by afmBaselineCorrection()
# via the afmDispatch() helper.
afmBaselineCorrectionSingle <-
  function(afmdata,
           ZPointApp = NULL,
           ZPointRet = NULL,
           fitpause = c("approach","retract","none"),
           vsTime = FALSE,
           sinusoidal = FALSE) {
    Segment <- Z <- NULL
    fitpause <- match.arg(fitpause)

    # If afmdata is a multi-indentation experiment, the contact point is
    # obtained with its first approach segment.
    approach_string <- ifelse(is.afmmulti(afmdata), "approach1", "approach")
    data <- afmdata$data
    isApproach <- data$Segment == approach_string

    if (is.null(ZPointApp)) {
      if ("CP" %in% names(afmdata)) {
        Zapp <- data$Z[isApproach]
        ZPointApp <- 0.3 * max(Zapp) + 0.7 * afmdata$CP[["CP"]]
      } else {
        stop(
          'ZPointApp should be given or, otherwise afmContactPoint() function
          should have been run first.'
        )
      }
    }

    # Baseline of the approach segment: best line fit on its off-contact
    # region (Z > ZPointApp).
    data.approach <- data[isApproach & data$Z > ZPointApp,
                          names(data) != "Segment", drop = FALSE]

    if (vsTime) {
      # Force vs time baseline (Moreno Flores, 2016): the fit is
      # subtracted from every segment of the curve.
      if (!sinusoidal) {
        fit.approach.time <- lm(Force ~ Time, data = data.approach)
      } else {
        fit.approach.time <- nls(Force ~ a0 + a1*Time + b * cos(w*Time+phi),
                                 data = data.approach,
                                 start = list(a0 = -0.5, a1 = -0.001, b = 0.001,
                                              w = 0.2/(2*pi), phi = pi/4))
      }
      afmdata$data$ForceCorrected <-
        data$Force -
        predict(fit.approach.time, data.frame(Time = data$Time))
      return(afmdata(afmdata))
    }

    # Force vs Z baseline: the approach and retract segments are corrected
    # with their own off-contact fits. The correction is assigned rowwise,
    # so the row order of the curve is preserved whatever its segment
    # layout, and every other segment (pause, modulation, ...) is
    # corrected according to the fitpause option.
    fit.approach <- lm(Force ~ Z, data = data.approach)
    ForceCorrected <- data$Force
    ForceCorrected[isApproach] <-
      data$Force[isApproach] -
      predict(fit.approach, data.frame(Z = data$Z[isApproach]))

    isRetract <- data$Segment == "retract"
    if (any(isRetract)) {
      Zret <- data$Z[isRetract]
      if (is.null(ZPointRet)) {
        if ("DP" %in% names(afmdata)) {
          ZPointRet <- 0.3 * max(Zret) + 0.7 * afmdata$DP[["DP"]]
        } else if ("CP" %in% names(afmdata)) {
          ZPointRet <- 0.6 * max(Zret) + 0.4 * afmdata$CP[["CP"]]
        } else{
          stop(
            'ZPointRet should be given or otherwise afmDetachPoint()
            and/or afmContactPoint() function(s) should have been run first.'
          )
        }
      }
      data.retract <- data[isRetract & data$Z > ZPointRet,
                           names(data) != "Segment", drop = FALSE]
      fit.retract <- lm(Force ~ Z, data = data.retract)
      ForceCorrected[isRetract] <-
        data$Force[isRetract] -
        predict(fit.retract, data.frame(Z = data$Z[isRetract]))
    }

    isMiddle <- !isApproach & !isRetract
    if (any(isMiddle) && fitpause != "none") {
      if (fitpause == "retract" && !any(isRetract)) {
        stop("fitpause = 'retract' requires a retract segment in the curve; ",
             "use 'approach' or 'none' instead.",
             call. = FALSE)
      }
      fit.middle <- if (fitpause == "retract") fit.retract else fit.approach
      ForceCorrected[isMiddle] <-
        data$Force[isMiddle] -
        predict(fit.middle, data.frame(Z = data$Z[isMiddle]))
    }
    afmdata$data$ForceCorrected <- ForceCorrected
    return(afmdata(afmdata))
  }
