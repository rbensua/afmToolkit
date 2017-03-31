#' @title Afmexperiment check.
#'   
#' @description Checks wether an R object is an afmexperiment or not.
#' @usage is.afmexperiment(x)
#' @param x Any \bold{R} object.
#' @return Returns \code{TRUE} if its argument is an afmdata (that is, has "afmexperiment"
#'   amongst its classes) and \code{FALSE} otherwise.
#'   
#' @export


is.afmexperiment <- function(x){
  inherits(x,"afmexperiment")
}