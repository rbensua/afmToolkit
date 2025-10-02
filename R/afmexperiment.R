#' @title AFM experiment
#'   
#' @description This function creates an \code{afmexperiment} structure, which is as list 
#' (or an array) of elements of \code{afmdata} class.
#' @usage afmexperiment(data, ID=NULL)
#' @param data A variable of  \code{afmdata} class, or a list of elements of
#'   \code{afmdata} class.
#' @param ID Character string with the identifier of the \code{data} variable or a string
#'   array in case \code{data} is a list of \code{afmdata} variables.
#'   
#' @return An object of class \code{afmexp}.
#'   
#' @examples
#' dataFolder <- paste(path.package("afmToolkit"), "afmexperiment",sep = "/")
#' dataFiles <- list.files(dataFolder, pattern = "force", full.names = FALSE)
#' data <- lapply(dataFiles, afmReadJPK, path = dataFolder)
#' names(data) <- dataFiles
#' data <- afmexperiment(data)
#' plot(data[[1]])
#' @seealso \code{\link{afmdata}}
#' @export
afmexperiment <-
  function(data, ID=NULL) {
    if (!is.afmexperiment(data)) {
        if (is.afmdata(data)){
          experiment <- list(data)
          if (is.null(ID)){
           names(experiment) <- "1" 
          }else {
          names(experiment) <- ID
          }
          return(structure(experiment, class = "afmexperiment"))
        }else if (is.list(data)){
          if (prod(sapply(data, is.afmdata))){
            experiment <- data
            if(is.null(ID)){
              if (is.null(names(data))){
                names(experiment) <- as.character(seq_along(data))
              }
            }else{
              names(experiment) <- ID
            }
          }else{
            stop("All elements in the list should be of afmdata class.")
          }
          return(structure(experiment, class = "afmexperiment"))
        }      
    } else{
      return(data)
    }
    
}