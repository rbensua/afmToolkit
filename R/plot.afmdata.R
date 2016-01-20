#' @title Plot an afmdata object
#' 
#' 
#' @description 
#' Plots an afmdata object.
#' 
#' @param \code{afmdata}: An object of \code{afmdata} class.
#' @param \code{vs}: The variable for the x-axis. May take the values "Time" or 
#'   "Z". It defaults to "Z", plotting thus a Force-Distance curve. If \code{vs}
#'   is set to "Time", then it plots  a Force-Time curve.
#' @param \code{segment}: The segment of the curve to be plotted. If
#'   \code{segment = "all"} then all segments of the curve are plotted. Possible
#'   values are: \code{"approach"}, \code{"contact"}, \code{"retract"} and 
#'   \code{"all"}.
#'   
#' 
#' @export
plot.afmdata <- function(afmdata, vs = "Z", segment = "all", ...)
{
  if (vs == "Time") {
    if (!("Time" %in% colnames(afmdata$data))) {
      stop("There is no column 'Time'")
    }
    if ("ForceCorrected" %in% names(afmdata$data)) {
      data <-
        subset(afmdata$data, select = c(Time,ForceCorrected,Segment))
      colnames(data) <- c("x","y","Segment")
    } else {
      data <- subset(afmdata$data, select = c(Time,Force,Segment))
      colnames(data) <- c("x","y","Segment")
    }
    if (segment == "all") {
      plt <- ggplot(data,aes(x = x,y = y, col = Segment)) + geom_line() +
        labs(x = "Time", y = "Force") + theme_bw()
    } else {
      data <- subset(data, Segment == segment)
      plt <- ggplot(data,aes(x = x,y = y)) + geom_line() +
        labs(x = "Time", y = "Force") + theme_bw()
    }
    if ("ForceCorrected" %in% names(afmdata$data) &
        segment != "contact") {
      plt +
        geom_hline(yintercept = 0,lty = 2)
    } else{
      plt
    }
  } else {
    if ("ForceCorrected" %in% names(afmdata$data)) {
      data <- subset(afmdata$data, select = c(Z,ForceCorrected,Segment))
      colnames(data) <- c("x","y","Segment")
    }else{
      data <- subset(afmdata$data, select = c(Z,Force,Segment))
      colnames(data) <- c("x","y","Segment")
    }
    if (segment == "all") {
      plt <- ggplot(data,aes(x = x,y = y, col = Segment)) + geom_line() +
        labs(x = "Distance", y = "Force") + theme_bw()
    }else{
      data <- subset(data, Segment == segment)
      plt <- ggplot(data,aes(x = x,y = y)) + geom_line() +
        labs(x = "Distance", y = "Force") + theme_bw()
    }
    if ("ForceCorrected" %in% names(afmdata$data) &
        segment != "contact") {
      plt +
        geom_hline(yintercept = 0,lty = 2)
    } else{
      plt
    }
  }
}