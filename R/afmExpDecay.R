#' @title Exponential decay fit
#'   
#'   
#' @description Fits a viscoelastic exponential decay in a Force-Relaxation or Creep
#'   experiments as described in Nanotechnology 2010 (see references).
#' @usage afmExpDecay(afmdata, nexp = 2, tmax = NULL, type = c("CH","CF"), plt = TRUE,
#'   ...)
#' @param afmdata An object of \code{afmdata} class with a \bold{pause} segment and a
#'   \bold{Time} column in the \code{data} dataframe.
#' @param nexp Number of expontials in the Prony series to be fitted. Currently only one
#'   or two exponentials are supported. Default is 2.
#' @param tmax Maximum time considered in the relaxation curve. It defaults to \code{Inf},
#'   meaning that the whole pause segment is considered.
#' @param type Type of the experiment. Can be either "CH" (Constant Height) for a
#'   force-relaxation experiment or "CF" (Constant Force) for a creep experiment. Default
#'   is \code{type = "CH"}.
#' @param plt Logical. If TRUE (default) then a plot of the pause segment with the overlay
#'   of the fit is shown.
#' @param ... Options passed to the \code{nlsM()} function from the \code{minpack.lm}
#'   package. At least should contain the starting values (\code{start = list(...)}) for
#'   the Levenberg-Mardquart nonlinear least square method.
#'   
#' @return An \code{afmdata} class variable which will consist on the original input
#'   \code{afmdata} variable plus a new list named \code{ExpFit} with the following
#'   fields:
#'   
#'   \code{expdecayModel}: A nls object returned from  \code{nlsM()} function.
#'   
#'   \code{expdecayFit}: The values predicted by the fit, returned from the 
#'   \code{predict()} function.
#' @examples 
#' data <- afmReadJPK("force-save-JPK-3h.txt.gz", path = path.package("afmToolkit"))
#' width <- 20
#' mul1 <- 1
#' mul2 <- 10
#' data <- afmContactPoint(data, width = width, mul1 = mul1, mul2 = mul2)
#' data <- afmDetachPoint(data, width = width, mul1 = mul1, mul2 = mul2)
#' data <- afmBaselineCorrection(data)
#' data <- afmExpDecay(data, nexp = 2, type = "CH")
#' @references Susana Moreno-Flores, Rafael Benitez, Maria dM Vivanco and Jose Luis 
#' Toca-Herrera (2010). "Stress relaxation and creep on living cells with the atomic force
#' microscope: a means to calculate elastic moduli and viscosities of cell components".
#' Nanotechnology, \strong{21} (44), pp. 445101.
#' @import minpack.lm
#' @import ggplot2
#' @importFrom stats predict
#' @export




afmExpDecay <- function(afmdata, nexp = 2, tmax = NULL,
                        type = c("CH","CF"), plt = TRUE, ...) {
  Segment <- Time <- Force <- Z <- NULL
  type <- match.arg(type)
  if (is.afmexperiment(afmdata)){
    afmdata <- lapply(afmdata, function(x){
      if (!is.null(x$params$curvename)) {
        print(paste("Processing curve: ", x$params$curvename), sep = " ")
      }
      afmExpDecay(x, nexp = nexp, tmax = tmax,type = type, plt = plt, ...)
      })
    return(afmexperiment(afmdata))
  }else if (!is.afmdata(afmdata)) {
    stop("Input must be an afmdata or an afmexperiment object!")
  }else{
  if (!"pause" %in% levels(afmdata$data$Segment)) {
    stop("A pause segment must be present!")
  }
  if (is.null(tmax)) {
    tmax <- Inf
  }
  if (type == "CF") {
    decay <-
      subset(afmdata$data, Segment == "pause" &
               Time <= tmax, select = c("Z", "Time"))
    decay$Time <- decay$Time - min(decay$Time)
    if (nexp == 1) {
      expdecayModel <-
        nlsLM(Z ~ c0 + c1 * exp(-x1 * Time) ,data = decay, ...)
    } else if (nexp == 2) {
      expdecayModel <-
        nlsLM(Z ~ c0 + c1 * exp(-x1 * Time) + c2 * exp(-x2 * Time),
              data = decay,...)
    } else {
      stop("No more than two exponentials are currently supported!")
    }
    expdecayFit <-
      predict(expdecayModel, data.frame(Time = decay$Time))
  }else if (type == "CH") {
    decay <- subset(
      afmdata$data, Segment == "pause" & Time <= tmax,
      select = c("ForceCorrected", "Time")
    )
    decay$Time <- decay$Time - min(decay$Time)
    if (nexp == 1) {
      expdecayModel <-
        nlsLM(ForceCorrected ~  a0 + a1 * exp(-Time / tau1),
              data = decay, ...)
    }else if (nexp == 2) {
      expdecayModel <-
        nlsLM(ForceCorrected ~ a0 + a1 * exp(-Time / tau1) +
                a2 * exp(-Time / tau2),
              data = decay,...)
    } else {
      stop("No more than two exponentials are currently supported!")
    }
    
    expdecayFit <-
      predict(expdecayModel, data.frame(Time = decay$Time))
  } else
    stop("type should be either 'CF' or 'CH'!")
  if (plt) {
    df <- subset(afmdata$data, Segment == "pause" & Time <= tmax,
                 select = c("Time"))
    if (type == "CH"){
    print(
      plot(afmdata, vs = "Time", segment = "pause") +
        geom_line(
          data = data.frame(Time = df, Force = expdecayFit),
          aes(x = Time, y = Force), col = "green", size = 1.5
        )
    )
    } else{
      pause <- subset(afmdata$data, Segment == "pause")
        print(ggplot(data = pause, aes(x = Time, y = Z)) + 
          geom_line() + 
          geom_line(
            data = data.frame(Time = df, Z = expdecayFit),
            aes(x = Time, y = Z), col = "green", size = 1.5
          ))
    }
  }
  ExpFit<- list(expdecayModel = expdecayModel, expdecayFit = expdecayFit) 
  return(append.afmdata(afmdata, ExpFit))
  }
}