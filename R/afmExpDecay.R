#' @title Exponential decay fit
#'   
#'   
#' @description Fits a viscoelastic exponential decay in a Force-Relaxation or
#' Creep experiments.
#' 
#' @param \code{afmdata}: An object of \code{afmdata} class with a 
#'   \bold{Contact} segment and a \bold{Time} column in the \code{data}
#'   dataframe.
#' @param \code{nexp}: Number of expontials in the Prony series to be fitted.
#'   Currently only one or two exponentials are supported. Default is 2.
#' @param \code{tmax}: Maximum time considered in the relaxation curve. It
#'   defaults to \code{Inf}, meaning that the whole contact segment is
#'   considered.
#' @param \code{type}: Type of the experiment. Can be either "CH" (Constant
#'   Height) for a force-relaxation experiment or "CF" (Constant Force) for a
#'   creep experiment. Default is \code{type = "CH"}.
#' @param \code{plt:} Logical. If TRUE (default) then a plot of the Contact
#'   segment with the overlay of the fit is shown.
#' @param \code{...}: Options passed to the \code{nlsM()} function from the
#'   \code{minpack.lm} package. At least should contain the starting values
#'   (\code{start = list(...)}) for the Levenberg-Mardquart nonlinear least
#'   square method.
#'   
#' @return A list of:
#'   
#'   \code{expdecayModel}: A nls object returned from  \code{nlsM()} function.
#'   
#'   \code{expdecayFit}: The values predicted by the fit, returned from the
#'   \code{predict()} function.
#'   
#' @export




afmExpDecay <- function(afmdata, nexp = 2, tmax = NULL,
                        type = c("CH","CF"), plt = TRUE, ...) {
  
  type <- match.arg(type)
  if (is.afmexperiment(afmdata)){
    expfit <- lapply(afmdata, function(x) afmExpDecay(x, nexp = nexp, tmax = tmax,
                                                      type = type, plt = plt, ...))
    afmdata <- mapply(afmdata,expfit, FUN = function(x,y) 
      append.afmdata(x,y,name = "ExpFit"), SIMPLIFY = FALSE)
    return(afmexperiment(afmdata))
  }else if (!is.afmdata(afmdata)) {
    stop("Input must be an afmdata or an afmexperiment object!")
  }else{
  if (!"contact" %in% levels(afmdata$data$Segment)) {
    stop("A contact segment must be present!")
  }
  if (is.null(tmax)) {
    tmax <- Inf
  }
  if (type == "CF") {
    decay <-
      subset(afmdata$data, Segment == "contact" &
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
      afmdata$data, Segment == "contact" & Time <= tmax,
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
    stop("type shoulf either 'CF' or 'CH'!")
  if (plt) {
    df <- subset(afmdata$data, Segment == "contact" & Time <= tmax,
                 select = c("Time"))
    if (type == "CH"){
    print(
      plot(afmdata, vs = "Time", segment = "contact") +
        geom_line(
          data = data.frame(Time = df, Force = expdecayFit),
          aes(x = Time, y = Force), col = "green", size = 1.5
        )
    )
    } else{
      contact <- subset(afmdata$data, Segment == "contact")
        print(ggplot(data = contact, aes(x = Time, y = Z)) + 
          geom_line() + 
          geom_line(
            data = data.frame(Time = df, Z = expdecayFit),
            aes(x = Time, y = Z), col = "green", size = 1.5
          ))
    }
  }
  return(list(expdecayModel = expdecayModel, expdecayFit = expdecayFit))
  }
}