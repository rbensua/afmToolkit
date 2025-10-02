#' @title Performs a baseline correction to an AFM F-z curve
#'
#' @description
#' This function performs the baseline correction to an AFM F-z curve within an
#' \code{afmdata} structure.
#'
#'  It substracts a best fit line to the cuve: for the approach and contact 
#'  segments, it fits a line to the approach curve
#' points where for which |z|>ZPointApp and for the retract segment,
#' it fits a line to the retract curve where |z|> ZpointRet.
#'
#' If no ZPointApp is given and the contact point has been already estimated 
#' (via \code{afmContactPoint()} function), then it is found as
#' \deqn{
#' ZPointApp = 0.7 ContactPoint + 0.3 max(Z)
#' }
#' @param afmdata An \code{afmdata} structure.
#' @param ZPointApp Point in the approach segment of the curve
#' that defines the approach baseline
#' @param ZPointRet Point in the retract segment of the curves that 
#' defines the retract baseline
#' @param fitpause Behaviour for the baseline correction at the pause segment: if "approach" (default), 
#' the pause segment is correted using the best line fit done on the approach segment, 
#' if "retract" the best line fit of the retract segment is used, if "none", no baseline correction 
#' is done on the pause segment. 
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
    Segment <- Z <- NULL
    fitpause <- match.arg(fitpause)
    if (is.afmexperiment(afmdata)){
      data <- lapply(afmdata, function(x) afmBaselineCorrection(x,
                                                                ZPointApp = ZPointApp,
                                                                ZPointRet = ZPointRet))
      return(afmexperiment(data))
    }else if (is.afmdata(afmdata)){
    
    # First determine how many segments are in the curve
    # If afmdata is a multi-indentation experiment, we obtain the contact point 
    # with the first approach segment.
      approach_string <- ifelse(is.afmmulti(afmdata),"approach1","approach")
      
      
    
    N <- nlevels(afmdata$data$Segment)
    if (is.null(ZPointApp)) {
      if ("CP" %in% names(afmdata)) {
        Zapp <- subset(afmdata$data, Segment == approach_string)$Z
        ZPointApp <- 0.3 * max(Zapp) + 0.7 * afmdata$CP[["CP"]]
      } else {
        stop(
          'ZPointApp should be given or, otherwise afmContactPoint() function
          should have been run first.'
        )
      }
    }
    
    data.approach <- subset(afmdata$data,
                            Segment == approach_string &
                              Z > ZPointApp,
                            select = -Segment)
    
    if(vsTime){
      if (!sinusoidal){
        fit.approach.time <- lm(Force  ~ Time, data = data.approach)
        afmdata$data$ForceCorrected <-
          afmdata$data$Force -
          predict(fit.approach.time, data.frame(Time = afmdata$data$Time))
      }else{
        print("Doing sinus")
        fit.approach.time <-nls(Force ~ a0 + a1*Time + b * cos(w*Time+phi), 
                                data = data.approach, 
                                start = list(a0 = -0.5, a1 = -0.001, b = 0.001, w = 0.2/(2*pi), phi = pi/4))
        afmdata$data$ForceCorrected <-
          afmdata$data$Force -
          predict(fit.approach.time, data.frame(Time = afmdata$data$Time))
      }
    }else{
      fit.approach <- lm(Force ~ Z, data = data.approach)  
      F.corrected.approach <-
        subset(afmdata$data, Segment == approach_string)$Force -
        predict(fit.approach, data.frame(Z = subset(afmdata$data,
                                                    Segment == approach_string)$Z))
      
      if (N == 1) {
        # If N = 1 there is only the approach segment.
        afmdata$data$ForceCorrected <- c(F.corrected.approach)
      } else if (N > 1) {
        # If N = 2 there are approach and retract segments.
        Zret <- subset(afmdata$data, Segment == "retract")$Z
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
        data.retract <- subset(afmdata$data,
                               Segment == "retract" &
                                 Z > ZPointRet,
                               select = -Segment)
        fit.retract <- lm(Force ~ Z, data = data.retract)
        F.corrected.retract <-
          subset(afmdata$data, Segment == "retract")$Force -
          predict(fit.retract, data.frame(Z = subset(afmdata$data,
                                                     Segment == "retract")$Z))
        if (N == 2) {
          afmdata$data$ForceCorrected <- c(F.corrected.approach,
                                           F.corrected.retract)
        } else{
         
          if (fitpause == "approach"){
            F.corrected.pause <-
              subset(afmdata$data, !Segment %in% c(approach_string,"retract"))$Force -
              predict(fit.approach, data.frame(Z = subset(afmdata$data,
                                                          !Segment %in% c(approach_string,"retract"))$Z))
          }else if (fitpause == "retract"){
            F.corrected.pause <-
              subset(afmdata$data, !Segment %in% c(approach_string,"retract"))$Force -
              predict(fit.retract, data.frame(Z = subset(afmdata$data,
                                                         Segment == "pause")$Z))
          }else{
            F.corrected.pause <- 
              subset(afmdata$data, Segment %in% c(approach_string,"retract"))$Force
          }
          afmdata$data$ForceCorrected <- c(F.corrected.approach,
                                           F.corrected.pause,
                                           F.corrected.retract)
        }
      }
    }
    return(afmdata(afmdata))
      }else{
    stop("input is not a valid afmdata or afmexperiment.")
      }
  }