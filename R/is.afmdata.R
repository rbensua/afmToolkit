#' @title Afmdata check.
#'   
#' @description Checks whether an R object is an afmdata or not.
#' @usage is.afmdata(x)
#' @param x Any \bold{R} object.
#' @return Returns \code{TRUE} if its argument is an afmdata (that is, has "afmdata" 
#'   amongst its classes) and \code{FALSE} otherwise.
#'   
#' @export


is.afmdata <- function(x){
  inherits(x,"afmdata")
}