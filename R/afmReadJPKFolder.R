#' @title Read all Nanowizard JPK ascii files in a folder
#'
#' @description
#' Read all JPK ascii files in a given folder. It searches for all files containing a 
#' given patter (".txt" by default) and uses the \code{afmReadJPJ} function.
#'
#' @usage afmReadJPKFolder(folder, pattern = ".txt", ...)
#' @param folder Name of the  folder containing the jpk files.
#' @param pattern Pattern that will identify the jok files (".txt" by default).
#' @param ... Other parameters passed to afmReadJPK function.
#' @return An \code{afmexperiment} class data structure with all F-d curves. 
#'
#' @examples
#' folder <- paste(path.package("afmToolkit"), "afmexperiment",sep = "/")
#' data <- afmReadJPKFolder(folder = folder)
#' str(data)
#' @export
#' 

afmReadJPKFolder <- function(folder, pattern = ".txt", ...){
  if (dir.exists(folder)){
  dataFiles <- list.files(folder, pattern = pattern, full.names = FALSE)
  data <- lapply(dataFiles, afmReadJPK, path = folder, ...)
  names(data) <- dataFiles
  return(afmexperiment(data))
  }else{
    stop("Folder does not exist...")
  }
}