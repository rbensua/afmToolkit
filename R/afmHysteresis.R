#' @title Hysteresis
#'
#' @description
#' Bla bla bal
#' 
#' @usage afmHysteresis(afmdata, silent = FALSE,...)
#' @param afmdata A Force-Distance curve with the afmdata structure. It should be a list
#'   with at least the 'data' field with a data frame of at least 4 columns.
#' @param silent Logical value. If TRUE omits the name of curve being processed (defaults to FALSE)
#' @param ... Extra options to be passed to \code{afmZeroPointSlope} function
#' @return An \code{afmdata} class variable which will consist on the original input
#'   \code{afmdata} variable plus a new list named \code{Hysteresis} with the following fields:
#'   
#' \code{Hyst_app} Time hyst app
#' 
#' \code{Hyst_ret} Time hyst ret
#'
#' \code{Hyst_tot} Time hyst total 
#'
#' \code{Hyst_ratio} Time hyst app/ret ratio
#' @examples
#' path <- path.package("afmToolkit")
#' data <- afmReadJPK("force-save-JPK-3h.txt.gz", path = path)
#' width <- 20
#' mul1 <- 1
#' mul2 <- 10
#' data <- afmContactPoint(data, width = width, mul1 = mul1, mul2 = mul2)
#' \dontrun{
#' plot(data, segment = "approach") + geom_vline(xintercept = data$CP$CP, lty = 2)
#' }
#' @seealso \code{\link{afmZeroPointSlope}}
#' @export



afmHysteresis <-
  function(afmdata, silent = FALSE,...) {
    Segment <- Time <-  NULL # For CRAN compatibility
    if(is.afmexperiment(afmdata)) {
      afmexperiment <- lapply(afmdata, function(x) {
        if (!is.null(x$params$curvename) & !silent) {
          print(paste("Processing curve: ", x$params$curvename), sep = " ")
        }
        afmHysteresis(x)
      })
      return(afmexperiment(afmexperiment))
    }else if(is.afmdata(afmdata)){
      if(!"Time" %in% colnames(afmdata$data)) {
        stop("A time column should be present!")
      }
      tmp1 <- afmZeroPointSlope(afmdata, segment = "approach", ...)
      tmp2 <- afmZeroPointSlope(afmdata, segment = "retract", ...)
      t1 <- tmp1$Slopes$t0Point
      t4 <- tmp2$Slopes$t0Point
      t2 <- max(afmdata$data$Time[afmdata$data$Segment == "approach"])
      t3 <- min(afmdata$data$Time[afmdata$data$Segment == "retract"])
      info = c(t1,t2,t3,t4)
      names(info) <- c("t1","t2","t3","t4")
      rm(list = c("tmp1","tmp2"))
      data.approach <- subset(afmdata$data, Segment == "approach" & Time >= t1)
      data.retract <- subset(afmdata$data, Segment == "retract" & Time <= t4)
      Hyst_app <- trapz(data.approach$Time, data.approach$ForceCorrected)
      Hyst_ret <- trapz(data.retract$Time, data.retract$ForceCorrected)
      Hyst_tot <- Hyst_app + Hyst_ret
      
      Energy_app <- abs(afmdata$params$speeds$approach) * Hyst_app
      Energy_ret <- abs(afmdata$params$speeds$retract) * Hyst_ret      
      Hysteresis <- Energy_app - Energy_ret
      Hyst_ratio <- Hysteresis / Energy_app
      Hysteresis <- list(
        Hyst_ratio  = Hyst_ratio,
        Hysteresis = Hysteresis,
        info = info
      )
      return(append.afmdata(afmdata,Hysteresis))

    }else{
      stop("No afmdata class input provided.")
    }
  }