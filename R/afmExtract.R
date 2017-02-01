#' @title Extract computed parameters from an afmexperiment
#'
#' @description
#' Extracts some parameters from an afmexperiment for an easy further analysis.
#' @usage afmExtract(afmexperiment, params = list("YM", "AE", "ED"))
#' @param afmexperiment Data of afmexperiment class.
#' @param params List of parameters to extract from the data.
#'
#' @return A data frame with the name of the curve and the corresponding values of the parameters extacted.
#'
#' @examples
#' data <- afmReadJPK("force-save-JPK-3h.txt",path = path.package("afmToolkit"))
#' str(data)
#' @export
#' 
afmExtract <- function(afmexperiment, params = list("YM", "AE", "ED")){
  if (is.afmexperiment(afmexperiment)){
   stop("Data should be of afmexperiment class!") 
  }
  extractedData <- data.frame(cuve = names(afmexperiment))
  if ("YM" %in% params){
    lapply(afmexperiment, function(x){ YM <- get("YoungModulus",get("YoungModulus",x))
    return(data.frame(YM = YM))})
    YM <- as.data.frame(do.call(rbind, YM), rownames = NULL)
    extractedData <- cbind(extractedData, YM)
  }
  if ("AE" %in% params){
    AE <- lapply(data, function(x){AE <- get("Energies", get("AdhEner",x))})
    AE <- as.data.frame(do.call(rbind, AE), rownames = NULL)
    extractedData <- cbind(extractedData, AE)
  }
  if ("ED" %in% params){
    stop("Not yet!")
  }
  return(extractedData)
}