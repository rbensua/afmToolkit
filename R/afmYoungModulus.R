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


afmYoungModulus <-
  function(afmdata,
           thickness = NULL,
           model = "Hertz",
           geometry = "pyramid",
           silent = TRUE,
           params) {
    if (is.afmexperiment(afmdata)) {
      afmexperiment <-
        lapply(afmdata, function(x){
          if(!is.null(x$params$curvename)){
            print(paste("Processing curve: ",x$params$curvename), sep = " ")
          }
          afmYoungModulus(
            x,
            thickness = thickness,
            model = model,
            geometry = geometry,
            params = params,
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
                 Segment == "approach" & Indentation < -6e-8 & Indentation > -8e-8)
      } else{
        fitdata <-
          subset(
            afmdata$data,
            Segment == "approach" & Indentation < 0 &
              abs(Indentation) < thickness
          )
      }
      #  fitdata$Indentation <- -fitdata$Indentation
      fitdata$Indentation <-
        fitdata$Indentation - fitdata$Indentation[1]
      fitdata$ForceCorrected <-
        fitdata$ForceCorrected - fitdata$ForceCorrected[1]
      fitLM <- lm(ForceCorrected  ~ Z, data = fitdata)
      fitYM <- lm(ForceCorrected ~ I(Indentation ^ 2) - 1, data = fitdata)
      fitYMfullPoly <-
        lm(ForceCorrected ~ poly(Indentation, 2, raw = TRUE), data = fitdata)
      predYM <-  predict(fitYM, newdata =
                           data.frame(Indentation = fitdata$Indentation))
      predYMFP <- predict(fitYMfullPoly,
                          newdata =
                            data.frame(Indentation = fitdata$Indentation))
      
      # plot(fitdata$Indentation,
      #      fitdata$ForceCorrected,
      #      xlab = "Indentation",
      #      ylab = "Force")
      # lines(fitdata$Indentation,
      #       predYM,
      #       col = "green",
      #       lwd = 2)
      # lines(fitdata$Indentation,
      #       predYMFP,
      #       col = "red",
      #       lwd = 2)
      # legend(
      #   "topright",
      #   c("Hertz - Sneddon", "Hertz+Adhesion"),
      #   lty = c(1, 1),
      #   col = c("green", "red")
      # )
      
      slope <- coef(fitYM)
      if (!silent) {
        print(summary(fitYM))
        print(summary(fitYMfullPoly))
      }
      if (is.null(params$nu)) {
        params$nu <- 0.5
      }
      YM <-
        as.numeric(slope * sqrt(2) * (1 - params$nu ^ 2) / tan(params$alpha * pi /
                                                                 180))
      #afmdata$params$YoungModulus <- YoungModulus
      #return(afmdata(afmdata))
      YoungModulus <- list(YoungModulus = YM,
                           fitYM = fitYM,
                           fitdata = subset(fitdata, 
                                            select = c(Indentation,ForceCorrected))) 
      return(append.afmdata(afmdata,YoungModulus))
    } else {
      stop("Error: input is not a valid afmdata or afmexperiment.")
    }
  }