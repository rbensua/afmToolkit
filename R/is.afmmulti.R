#' @title Afmdata multiindentation check.
#'   
#' @description Checks whether an R object is of class afmmulti or not.
#' @usage is.afmmulti(x)
#' @param x Any \bold{R} object.
#' @return Returns \code{TRUE} if its argument is an afmmulti (that is, has "afmmulti" 
#'   amongst its classes) and \code{FALSE} otherwise.
#'   
#' @export


is.afmmulti <- function(x){
  inherits(x,"afmmulti")
}