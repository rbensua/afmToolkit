#' @title Plot an afmdata object
#'   
#'   
#' @description Plots an afmdata object.
#' 
#' @param x An object of \code{afmdata} class.
#' @param y Variable added for compatibility with \code{plot}.
#' @param vs The variable for the x-axis. May take the values "Time" or "Z". It defaults
#'   to "Z", plotting thus a Force-Distance curve. If \code{vs} is set to "Time", then it
#'   plots  a Force-Time curve.
#' @param segment The segment of the curve to be plotted. If \code{segment = "all"} then
#'   all segments of the curve are plotted. Possible values are: \code{"approach"},
#'   \code{"pause"}, \code{"retract"} and \code{"all"}.
#' @param ... Additional parameters to be pased to the ggplot functions.
#'   
#' @examples
#' # Loading the data
#' path <- path.package("afmToolkit")
#' data <- afmReadJPK("force-save-JPK-3h.txt.gz", path = path)
#' # Standard plot (out of the box)
#' plot(data)
#' # Computing the contact and detach points
#' data <- afmContactPoint(data, width = 20, mul1 = 1, mul2 = 10)
#' data <- afmDetachPoint(data, width = 40, mul1 = 3, mul2 = 20)
#' # Making the baseline correction
#' data <- afmBaselineCorrection(data)
#' # Plot once the baseline correction is done
#' plot(data)
#' # Plotting only retract segment
#' plot(data, segment = "retract")
#' # Plotting the pause segment: Force vs Time
#' plot(data, segment = "pause", vs = "Time")
#' @method plot afmdata
#' @export
plot.afmdata <- function(x, y = NULL, vs = "Z", segment = "all", ...)
{
  Segment <- Time <- ForceCorrected <- Force <- Z <- NULL
  if (!is.null(y)){
    y <- NULL
  }
  afmdata <- x
  if (is.afmdata(afmdata)){
    segment <- match.arg(segment, c("approach","pause","retract","all"))
   # vs <- match.arg(vs, c("Z","Time"))
  if (vs == "Time") {
    if (!("Time" %in% colnames(afmdata$data))) {
      stop("There is no column 'Time'")
    }
    if ("ForceCorrected" %in% names(afmdata$data)) {
      data <-
        subset(afmdata$data, select = c(vs,"ForceCorrected","Segment"))
      colnames(data) <- c("x","y","Segment")
    } else {
      data <- subset(afmdata$data, select = c(vs,"Force","Segment"))
      colnames(data) <- c("x","y","Segment")
    }
    if (segment == "all") {
      plt <- ggplot(data,aes(x = x,y = y, col = Segment)) + geom_line() +
        labs(x = "Time", y = "Force") + ggtitle(afmdata$params$curvename) + theme_bw()
    } else {
      data <- subset(data, Segment == segment)
      plt <- ggplot(data,aes(x = x,y = y)) + geom_line() +
        labs(x = "Time", y = "Force") + ggtitle(afmdata$params$curvename) +
        theme_bw()
    }
    if ("ForceCorrected" %in% names(afmdata$data) &
        segment != "pause") {
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
        labs(x = "Distance", y = "Force")  + ggtitle(afmdata$params$curvename) +
        theme_bw()
    }else{
      data <- subset(data, Segment == segment)
      plt <- ggplot(data,aes(x = x,y = y)) + geom_line() +
        labs(x = "Distance", y = "Force")  + ggtitle(afmdata$params$curvename) +
        theme_bw()
    }
    if ("ForceCorrected" %in% names(afmdata$data) &
        segment != "pause") {
      plt +
        geom_hline(yintercept = 0,lty = 2)
    } else{
      plt
    }
  }
  }else{
    stop("Input is not an afmdata class variable!")
  }
}