#' @title Stress relaxation/ Creep model fit
#'   
#'   
#' @description Fits a viscoelastic exponential decay in a Force-Relaxation or Creep
#'   experiments as described in Nanotechnology 2010 (see references).
#' @usage afmRelax(afmdata, model = c("power","exp"), nexp = 2, tmax = NULL, 
#'        type = c("CH","CF"), plt = TRUE, ...)
#' @param afmdata An object of \code{afmdata} class with a \bold{pause} segment and a
#'   \bold{Time} column in the \code{data} dataframe.
#' @param model Type of model to be fitted. Could be either "power" for a power 
#' law model or "exp" for an exponential decay model. 
#' @param nexp Number of expontials in the Prony series to be fitted if \code{model} is \code{"exp"}. 
#' Currently only one or two exponentials are supported. Default is 2.
#' @param tmax Maximum time considered in the relaxation curve. It defaults to \code{Inf},
#'   meaning that the whole pause segment is considered.
#' @param type Type of the experiment. Can be either "CH" (Constant Height) for a
#'   force-relaxation experiment or "CF" (Constant Force) for a creep experiment. Default
#'   is \code{type = "CH"}.
#' @param plt Logical. If TRUE (default) then a plot of the pause segment with the overlay
#'   of the fit is shown.
#' @param ... Options passed to the \code{nlsM()} function from the \code{minpack.lm}
#'   package. At least should contain the starting values (\code{start = list(...)}) for
#'   the Levenberg-Mardquart nonlinear least square method. For use only if \code{model} is \code{"exp"}. 
#'   
#' @return An \code{afmdata} class variable which will consist on the original input
#'   \code{afmdata} variable plus a new list named \code{RelaxFit} with the following
#'   fields:
#'   
#'   \code{model}: The model used in the fit.
#'   
#'   \code{relaxModel}: A nls object returned from  \code{nlsM()} function.
#'   
#'   \code{relaxFit}: The values predicted by the fit, returned from the 
#'   \code{predict()} function.
#' @examples 
#' data <- afmReadJPK("force-save-JPK-3h.txt.gz", path = path.package("afmToolkit"))
#' width <- 20
#' mul1 <- 1
#' mul2 <- 10
#' data <- afmContactPoint(data, width = width, mul1 = mul1, mul2 = mul2)
#' data <- afmDetachPoint(data, width = width, mul1 = mul1, mul2 = mul2)
#' data <- afmBaselineCorrection(data)
#' data <- afmRelax(data, model = "exp", nexp = 2, type = "CH")
#' @references Susana Moreno-Flores, Rafael Benitez, Maria dM Vivanco and Jose Luis 
#' Toca-Herrera (2010). "Stress relaxation and creep on living cells with the atomic force
#' microscope: a means to calculate elastic moduli and viscosities of cell components".
#' Nanotechnology, \strong{21} (44), pp. 445101.
#' @import minpack.lm
#' @import ggplot2
#' @importFrom stats predict
#' @export




afmRelax <- function(afmdata, model = c("power","exp"), nexp = 2, tmax = NULL,
                     type = c("CH","CF"), plt = TRUE, ...) {
  Segment <- Time <- Z <- Force <- ForceCorrected <- Force_fitted <- 
    Z_fitted <- NULL # For CRAN compatibility
  model <- match.arg(model)
  type <- match.arg(type)
  if (is.afmexperiment(afmdata)){
    afmdata <- lapply(afmdata, function(x){
      if (!is.null(x$params$curvename)) {
        print(paste("Processing curve: ", x$params$curvename), sep = " ")
      }
      afmRelax(x, model = model, nexp = nexp, tmax = tmax,type = type, plt = plt, ...)
    })
    return(afmexperiment(afmdata))
  }else if (!is.afmdata(afmdata)) {
    stop("Input must be an afmdata or an afmexperiment object!")
  }else{
    if (!any(grepl("pause", levels(afmdata$data$Segment)))) {
      stop("A pause segment must be present!")
    }
    if (is.null(tmax)) {
      tmax <- Inf
    }
    if (!is.afmmulti(afmdata)){
      seg <- "pause"
      fitting <- fitcurve(afmdata, model = model, type = type, segment = seg, tmax = tmax, nexp = nexp,...)
      Relax <- list(model = model, relaxModel = fitting$relaxModel, relaxFit = fitting$relaxFit) 
    } else{
      seg <- as.character(unique(afmdata$data$Segment)[grepl("pause",unique(afmdata$data$Segment))])
      fitting <- lapply(seg, function(x) fitcurve(afmdata, model = model, type = type, segment = x, tmax = tmax, nexp = nexp,...))
      names(fitting) <- seg  
      Relax <- list(model = model, 
                    relaxModel = lapply(fitting, function(x) x$relaxModel), 
                    relaxFit = lapply(fitting, function(x) x$relaxFit)) 
    }
    
    
    if (plt) {
      if(!is.afmmulti(afmdata)){
        df <- subset(afmdata$data, Segment == "pause" & Time <= tmax,
                     select = c("Time"))
        if (type == "CH"){
          print(
            #plot(afmdata, vs = "Time", segment = "pause") +
              ggplot(data = subset(afmdata$data, Segment == "pause"), 
                     aes(x = Time, y = ForceCorrected)) +
              geom_point(alpha = 0.5, shape = 21, colour = "black", fill = "lightblue") + 
              geom_line(
                data = data.frame(Time = df, Force = Relax$relaxFit),
                aes(x = Time, y = Force), col = "blue", lwd = 1) +
              theme_bw() + ylab("Force (N)") + ggtitle(afmdata$params$curvename)
          )
        } else{
          pause <- subset(afmdata$data, Segment == "pause")
          print(ggplot(data = pause, aes(x = Time, y = Z)) + 
                  geom_point(alpha = 0.4, shape = 21, colour = "black", fill = "lightblue") + 
                  geom_line(
                    data = data.frame(Time = df, Z = Relax$relaxFit),
                    aes(x = Time, y = Z), col = "blue", lwd = 1) + 
                  theme_bw() + ylab("Z (m)") + ggtitle(afmdata$params$curvename)
                )
        }
      }else{
        df <- afmdata$data %>% filter(grepl("pause", Segment) & Time <= tmax)
        if (type == "CH"){
          df <- df %>% mutate(Force_fitted = unlist(Relax$relaxFit))
          print(ggplot(data = df, aes(x = Time, y = ForceCorrected, fill = Segment)) + 
                  geom_point(alpha = .5, shape = 21, colour = "black") +
                  geom_line(aes(y = Force_fitted, group = Segment, col = Segment), lwd = 1.1) +
                  theme_bw() + ylab("Force (N)") + ggtitle(afmdata$params$curvename))  
            
        } else{
          df <- df %>% mutate(Z_fitted = unlist(Relax$relaxFit))
          print(ggplot(data = df, aes(x = Time, y = Z, fill = Segment)) + 
                  geom_point(alpha = .5, shape = 21, colour = "black") +
                  geom_line(aes(y = Z_fitted, group = Segment, col = Segment), lwd = 1.1) +
                  theme_bw() + ylab("Z (m)") + ggtitle(afmdata$params$curvename))
        }
      }
    } 
    return(append.afmdata(afmdata, Relax))
  }
}


fitcurve <- function(afmdata, type, model, segment, tmax, nexp = 2,...){
  Segment <- Time <- NULL
  if (type == "CF") {
    decay <-
      subset(afmdata$data, Segment == segment &
               Time <= tmax, select = c("Z", "Time"))
    decay$Time <- decay$Time - min(decay$Time) + 0.1*diff(decay$Time)[1]
    if (model == "exp"){
      if (nexp == 1) {
        relaxModel <-
          nlsLM(Z ~ c0 + c1 * exp(-x1 * Time) ,data = decay, ...)
      } else if (nexp == 2) {
        relaxModel <-
          nlsLM(Z ~ c0 + c1 * exp(-x1 * Time) + c2 * exp(-x2 * Time),
                data = decay,...)
      } else {
        stop("No more than two exponentials are currently supported!")
      }
      relaxFit <-
        predict(relaxModel, data.frame(Time = decay$Time))
    } else {
      relaxModel <- lm( log(Z) ~ log(Time), data = decay )
      relaxFit <-
        exp(predict(relaxModel, data.frame(Time = decay$Time)))
    }
    
    
  }else if (type == "CH") {
    
    decay <- subset(
      afmdata$data, Segment == segment & Time <= tmax,
      select = c("ForceCorrected", "Time")
    )
    
    if (model == "exp"){
      decay$Time <- decay$Time - min(decay$Time) 
      if (nexp == 1) {
        relaxModel <-
          nlsLM(ForceCorrected ~  a0 + a1 * exp(-Time / tau1),
                data = decay, ...)
      }else if (nexp == 2) {
        relaxModel <-
          nlsLM(ForceCorrected ~ a0 + a1 * exp(-Time / tau1) +
                  a2 * exp(-Time / tau2),
                data = decay,...)
      } else {
        stop("No more than two exponentials are currently supported!")
      }
      relaxFit <-
        predict(relaxModel, data.frame(Time = decay$Time))
    } else {
      relaxModel <- lm( log(ForceCorrected) ~ log(Time), data = decay )
      relaxFit <-
        exp(predict(relaxModel, data.frame(Time = decay$Time)))
    }
    
    
  }
  
  return(list(relaxModel = relaxModel, relaxFit = relaxFit))
}

