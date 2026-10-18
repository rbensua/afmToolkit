#' @title Dispatch a function over an afmdata or an afmexperiment
#'
#' @description
#' Internal helper used by the analysis functions of the package
#' (\code{afmContactPoint()}, \code{afmDetachPoint()}, \code{afmBaselineCorrection()},
#' \code{afmZeroPointSlope()}, \code{afmAdhesionEnergy()}, ...).
#'
#' Given an input that is either a single \code{afmdata} curve or an
#' \code{afmexperiment} (a list of curves), it applies \code{FUN} to the
#' single curve or to each curve of the experiment and returns the result
#' with the corresponding class. It also centralises the
#' "Processing curve" message and the error raised when the input is of
#' neither class.
#' @param afmdata An \code{afmdata} or \code{afmexperiment} variable.
#' @param FUN Internal (single-curve) function to be applied.
#' @param ... Additional arguments passed to \code{FUN}.
#' @param silent Logical. If TRUE, suppresses the "Processing curve" message
#' when \code{afmdata} is an \code{afmexperiment}. Defaults to FALSE.
#' @param error Error message raised when \code{afmdata} is neither an
#' \code{afmdata} nor an \code{afmexperiment}.
#' @return The result of \code{FUN(afmdata, ...)} for a single curve, or an
#' \code{afmexperiment} with \code{FUN} applied to each of its curves.
#' @keywords internal
afmDispatch <-
  function(afmdata,
           FUN,
           ...,
           silent = FALSE,
           error = "No afmdata or afmexperiment class input provided.") {
    if (is.afmexperiment(afmdata)) {
      n_curves <- length(afmdata)
      ID <- names(afmdata)
      afmexperiment <- lapply(seq_len(n_curves), function(idx) {
        x <- afmdata[[idx]]
        if (!is.null(x$params$curvename) & !silent) {
          cat("Processing curve", idx , "/", n_curves,":", x$params$curvename,"\n")
        }
        FUN(x, ...)
      })
      return(afmexperiment(afmexperiment, ID = ID))
    } else if (is.afmdata(afmdata)) {
      return(FUN(afmdata, ...))
    } else {
      stop(error)
    }
  }
