#' @title Read all Bruke Nanoscope Veeco ascii files in a folder
#'
#' @description
#' Read all Veeco ascii files in a given folder. It searches for all files containing a 
#' given patter (".txt" by default) and uses the \code{afmReadVeeco} function.
#'
#' @usage afmReadVeecoFolder(folder, pattern = ".txt", ...)
#' @param folder Name of the  folder containing the Veeco files.
#' @param pattern Pattern that will identify the Veeco files (".txt" by default).
#' @param ... Parameters to be passed to the afmReadVeeco() function.

#' @return An \code{afmexperiment} class data structure with all F-d curves. 
#'
#' @examples
#' folder <- paste(path.package("afmToolkit"), "veecoFolder",sep = "/")
#' data <- afmReadVeecoFolder(folder = folder)
#' str(data)
#' @export
#' 

afmReadVeecoFolder <- function(folder, pattern = ".txt",...){
  if (dir.exists(folder)){
    dataFiles <- list.files(folder, pattern = pattern, full.names = FALSE)
    data <- lapply(dataFiles, afmReadVeeco, path = folder,...)
    names(data) <- dataFiles
    return(afmexperiment(data))
  }else{
    stop("Folder does not exist...")
  }
}