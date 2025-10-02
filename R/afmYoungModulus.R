#' @title afmYoungModulus
#'   
#' @description This function computes the Young's Modulus of the sample from the approach
#'   curve using Hertz's  contact model  for a pyramidal, paraboloid or conical tip.
#' @param afmdata An \code{afmdata} object. It should be a valid afmdata object upon which
#'   the Contact Point, the baseline correction and the Zero Force Point and the 
#'   Indentation must have been calculated first (using functions 
#'   \code{afmContactPoint()}, \code{afmBaselineCorrection()}, \code{afmZeroPointSlope()},
#'   and \code{afmIndentation()})
#' @param thickness Thickness (in m) of the surface. The Force - Indentation fit will be 
#'   done for values of the Indentation variable smaller than the thickness. If no value 
#'   is given, it will be done for all values in the curve for which the Indentation is 
#'   negative.
#' @param model Contact mechanics model to be used. Currently only Hertz's pure elastic 
#'   model is available.
#' @param geometry Geometry of the tip. Currently only pyramidal (default), paraboloid, 
#'   and conical geometries are implemented.
#' @param thin_correction Logical value. If TRUE a thin layer correction is performed as 
#' shown in Dimitriadis et al. 2002, and in Gavara and Chadwick 2012. 
#' @param showfit Logical value. If FALSE it prints the fit model summary (via 
#'   \code{summary.lm()}). Default value is FALSE
#' @param params A list containing different parameters of the model: e.g. nu (Poisson's 
#'   ratio) or alpha (internal angle, in degrees, of the pyramidal and conical tips) or R (tip radius, 
#'   in the paraboloid geometry)
#' @param silent Logical value. If TRUE suppresses the message indicating the name of the curve being processed
#' (useful for batch-processing large number of curves). Defaults to FALSE
#' @return An \code{afmdata} class variable which will consist on the original input 
#'   \code{afmdata} variable plus a new list named \code{YoungModulus} with the following 
#'   fields:
#'   
#'   \code{YoungModulus} The Young's modulus value (in Pa).
#'   
#'   \code{fitYM} The Force vs Indentation^2 fit as an \code{lm} object.
#'   
#'   \code{fitdata} The subset of the data used in the fit.
#' @examples 
#' data <- afmReadJPK("force-save-JPK-2h.txt.gz", path = path.package("afmToolkit"))
#' data <- afmContactPoint(data, width = 20, mul1 = 1, mul2 = 20)
#' data <- afmDetachPoint(data, width = 40, mul1 = 3, mul2 = 40)
#' data <- afmBaselineCorrection(data)
#' data <- afmZeroPointSlope(data, segment = "approach")
#' data <- afmIndentation(data)
#' data <- afmYoungModulus(data, thickness = 1e-8, params = list(alpha = 22),
#'                         silent = TRUE)
#' print(data$YoungModulus$YoungModulus)
#' @importFrom stats coef lm predict
#' @export


afmYoungModulus <-
  function(afmdata,
           thickness = NULL,
           model = "Hertz",
           geometry = c("pyramid","paraboloid", "cone"),
           thin_correction = FALSE,
           showfit = FALSE,
           params, 
           silent = FALSE) {
    Segment <- Indentation <- ForceCorrected <- NULL
    if (is.afmexperiment(afmdata)) {
      afmexperiment <-
        lapply(afmdata, function(x){
          if(!is.null(x$params$curvename)& !silent){
            print(paste("Processing curve: ",x$params$curvename), sep = " ")
          }
          afmYoungModulus(
            x,
            thickness = thickness,
            model = model,
            geometry = geometry,
            params = params,
            thin_correction = thin_correction,
            showfit = showfit,
            silent = silent)
          })
      return(afmexperiment(afmexperiment))
    } else if (is.afmdata(afmdata)) {
      if (!("ForceCorrected" %in% names(afmdata$data))) {
        stop("Baseline correction should be done first!")
      }
      if (!any ("Slopes" %in% names(afmdata))) {
        stop("Zero Point Force should be found first (run afmSlopes)!")
      }
      if (!afmdata$params$SpringConstant > 0) {
        stop("The cantilever spring constant should be provided!")
      }
      if (!("Indentation" %in% names(afmdata$data))) {
        stop("Indentation should be computed first(run afmIndentation)!")
      }
      if (is.null(thickness)) {
        fitdata <-
          subset(afmdata$data,
                 Segment == "approach" & Indentation <= 0)
      } else{
        fitdata <-
          subset(
            afmdata$data,
            Segment == "approach" & Indentation <= 0 &
              abs(Indentation) < thickness
          )
      }
      #  fitdata$Indentation <- -fitdata$Indentation
      fitdata$Indentation <-
        fitdata$Indentation - fitdata$Indentation[1]
      fitdata$ForceCorrected <-
        fitdata$ForceCorrected - fitdata$ForceCorrected[1]
      geometry <- match.arg(geometry)
      if (geometry == "pyramid"){
        if ( nrow(fitdata) > 3 ){
          if(thin_correction){
            if(is.null(params$h)){
              stop("You need to provide the height of the sample!")
            }
            warning("No thin sample correction for pyramidal geometry is available. Switching to conical geometry!")
            fitdata$Indentation <- -fitdata$Indentation
            fitdata <- subset(fitdata, Indentation>=0)
            chi <- fitdata$Indentation / params$h
            fitdata$corr_factor <- 1 + 1.7795 * 2 * tan(params$alpha * pi / 180) * chi / pi^2 + 
              16 * 1.7795^2 * tan(params$alpha * pi / 180)^2 * chi^2 
            fitYM <- lm(ForceCorrected ~ I(corr_factor * Indentation ^ 2) - 1, data = fitdata)
          } else{
            fitYM <- lm(ForceCorrected ~ I(Indentation ^ 2) - 1, data = fitdata)
          }
          predYM <-  predict(fitYM, newdata =
                               data.frame(Indentation = fitdata$Indentation))
          
          slope <- coef(fitYM)
          if (showfit) {
            print(summary(fitYM))
          }
          
          if(thin_correction){
          YM <-
            as.numeric(slope * 3 * pi / (8 *  
                         tan(params$alpha * pi / 180))
                       )
          } else {
            if (is.null(params$nu)) {
              params$nu <- 0.5
            }
            YM <-
              as.numeric(slope * sqrt(2) * (1 - params$nu ^ 2) / 
                           tan(params$alpha * pi / 180))
          }
        } else {
          YM <- NA
          fitYM <- NA
        }
      }else if (geometry == "cone"){
        if ( nrow(fitdata) > 3 ){
          if (thin_correction){
            if(is.null(params$h)){
              stop("You need to provide the height of the sample!")
            }
            fitdata$Indentation <- -fitdata$Indentation
            fitdata <- subset(fitdata, Indentation>=0)
            chi <- fitdata$Indentation / params$h
            fitdata$corr_factor <- 1 + 1.7795 * 2 * tan(params$alpha * pi / 180) * chi / pi^2 + 
              16 * 1.7795^2 * tan(params$alpha * pi / 180)^2 * chi^2 
            
            fitYM <- lm(ForceCorrected ~ I(corr_factor * Indentation ^ 2) - 1, data = fitdata)
            predYM <-  predict(fitYM)
            
            slope <- coef(fitYM)
            YM <-
              as.numeric(slope * 3 * pi / (8 *  
                                             tan(params$alpha * pi / 180))
              )
          }else {
            fitYM <- lm(ForceCorrected ~ I(Indentation ^ 2) - 1, data = fitdata)
            
            predYM <-  predict(fitYM, newdata =
                                 data.frame(Indentation = fitdata$Indentation))
            
            slope <- coef(fitYM)
            if (showfit) {
              print(summary(fitYM))
            }
            if (is.null(params$nu)) {
              params$nu <- 0.5
            }
            YM <-
              as.numeric(slope * pi * (1 - params$nu ^ 2) /
                           (2* tan(params$alpha * pi /180)))
          }
        } else { 
          YM <- NA
          fitYM <- NA
        }
        
      }else {
        fitdata$Indentation <- -fitdata$Indentation
        fitdata <- subset(fitdata, Indentation>=0)
        if (nrow(fitdata) > 3){
          if(thin_correction){
            params$nu <- 0.5
            
            chi <- sqrt(params$R * fitdata$Indentation) / params$h
            fitdata$corr_factor <- 1 + 1.133 * chi + 1.283 * chi^2 + 
              0.769 * chi^3 + 0.0975 * chi^4
            fitYM <- lm(ForceCorrected ~ I(corr_factor * Indentation^1.5) - 1, data = fitdata)
            
            predYM <-  predict(fitYM)
            
            slope <- coef(fitYM)
            if (!silent) {
              print(summary(fitYM))
            }
            
            YM <-
              as.numeric(slope * 9 / (16 * sqrt(params$R)))
          } else{
            fitYM <- lm(ForceCorrected ~ I(Indentation ^ 1.5) - 1, data = fitdata)
            
            predYM <-  predict(fitYM)
            
            slope <- coef(fitYM)
            if (!silent) {
              print(summary(fitYM))
            }
            if (is.null(params$nu)) {
              params$nu <- 0.5
            }
            YM <-
              as.numeric(slope * 0.75 * (1 - params$nu ^ 2) / sqrt(params$R))
          }
        } else {
          YM <- NA
          fitYM <- NA
        }
      }
      YoungModulus <- list(YoungModulus = YM,
                           fitYM = fitYM,
                           fitdata = subset(fitdata, 
                                            select = c(Indentation,ForceCorrected))) 
      return(append.afmdata(afmdata,YoungModulus))
    } else {
      stop("Error: input is not a valid afmdata or afmexperiment.")
    }
  }