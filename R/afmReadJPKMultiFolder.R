#' @title Read all Nanowizard Multi-indentation JPK ascii files in a folder
#'
#' @description
#' Read all multi-indentation JPK ascii files in a given folder. It searches for all files containing a 
#' given patter (".txt" by default) and uses the \code{afmReadJPJ} function.
#'
#' @usage afmReadJPKMultiFolder(folder, pattern = ".txt", ...)
#' @param folder Name of the  folder containing the jpk files.
#' @param pattern Pattern that will identify the jok files (".txt" by default).
#' @param ... Other parameters passed to afmReadJPKMultiIndent function.
#' @return An \code{afmmultiexp} class data structure with all F-d curves. 
#'
#' @examples 
#' folder <- paste(path.package("afmToolkit"), "multiIndentExperiment",sep = "/")
#' data <- afmReadJPKMultiFolder(folder = folder)
#' 
#' @export

afmReadJPKMultiFolder <- function(folder, pattern = ".txt", ...){
  if (dir.exists(folder)){
  dataFiles <- list.files(folder, pattern = pattern, full.names = FALSE)
  data <- lapply(dataFiles, afmReadJPKMultiIndent, path = folder, ...)
  names(data) <- dataFiles
  return(structure(data, class = c("afmexperiment","afmmultiexp")))
  }else{
    stop("Folder does not exist...")
  }
}