#' @title Append to an \code{afmdata} list.
#'   
#' @description This function appends a list to an existing afmdata structure. It is used
#' internally by several afm* functions when attaching the results to the input afmdata
#' variable. This function should not be used directly unless by experienced users.
#' 
#' @param afmdata   The afmdata to which the new list is going to be joined.
#' @param x A list to be appended.
#' @param name The name of new field of the resulting afmdata object. If none is given, it
#'   is the same as \code{x}.
#' @return The new list of class \code{afmdata}
#' @export

append.afmdata <- function(afmdata,x, name = NULL){
  classes <- class(afmdata)
  if (is.null(name)) name = deparse(substitute(x))
  if (name %in% names(afmdata)){
    afmdata[[name]] <- x
  }else {
  newlist <- list()
  newlist[[name]] <- x
  afmdata <- append(afmdata, newlist)
  }
  res <- afmdata(afmdata)
  class(res) <- classes
  return(res)
}