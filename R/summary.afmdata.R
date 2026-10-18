#' @title Summary of an \code{afmdata} class object.
#'   
#'   
#' @description This function sumarizes the main features of an afmdata object and,
#' optionnaly plots all segments available with all parameters estimated.
#' 
#' @param object An object of \code{afmdata} class.
#' @param plt Logical variable. If TRUE plots all available segments with all available
#'   data.
#' @param ... Additional arguments (for compatibility with \code{summary})
#' @examples
#' \dontrun{path <- path.package("afmToolkit")
#' data <- afmReadJPK("force-save-JPK-3h.txt.gz", path = path)
#' data <- afmContactPoint(data, width = 20, mul1 = 1, mul2 = 10)
#' data <- afmDetachPoint(data, width = 20, mul1 = 2, mul2 = 30)
#' data <- afmBaselineCorrection(data)
#' data <- afmAdhesionEnergy(data, width = 20, mul = 10)
#' data <- afmZeroPointSlope(data, segment = "approach")
#' data <- afmIndentation(data)
#' data <- afmYoungModulus(data, thickness = 1e-7, params = list(alpha = 22),
#'                         silent = TRUE)
#' data <- afmRelax(data, model = "exp", nexp = 2, type = "CH")
#' summary(data)}                        
#' @import grid
#' @importFrom  gridExtra ttheme_default tableGrob
#' @import scales
#' @import dplyr
#' @importFrom graphics curve plot
#' @importFrom stats coefficients complete.cases
#' @importFrom utils read.table
#' @rdname summary
#' @export
#' @method summary afmdata

# The summary method is now a thin orchestrator over internal modules:
#   afmSummaryTable()        - numeric (console) summary
#   afmSummaryApproachPlot() - approach segment plot
#   afmSummaryPausePlot()    - pause segment plots (relaxation fit)
#   afmSummaryRetractPlot()  - retract segment plot (adhesion energies)
summary.afmdata <- function(object, plt = TRUE, ...){
  afmdata <- object
  if (!is.afmdata(afmdata)){
    stop("Input should be an afmdata class object!")
  }
  res <- afmSummaryTable(afmdata)
  if (plt) {
    afmSummaryApproachPlot(afmdata, res$output)
    afmSummaryPausePlot(afmdata, res$typeoffit)
    afmSummaryRetractPlot(afmdata)
  }
  invisible(res$output)
}
