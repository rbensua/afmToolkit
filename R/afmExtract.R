#' @title Extract computed parameters from an afmexperiment
#'   
#' @description Extracts some parameters from an afmexperiment for an easy further
#' analysis.
#' @usage afmExtract(afmexperiment, params = list("YM", "AE", "ED"), opt.param = NULL)
#' @param afmexperiment Data of afmexperiment class.
#' @param params List of parameters to extract from the data.
#' @param opt.param Optional parameter or factor in the params field of the afmdata list 
#'   to add to the data extraction.
#' @return A data frame with the name of the curve and the corresponding values of the 
#'   parameters extacted.
#'   
#' @examples
#' \dontrun{
#' require(dplyr) # Not really necessary
#' 
#' # Load the data
#' data(batchExperiment)
#' 
#' # Process the afmexperiment
#' data <- afmContactPoint(batchExperiment, width = 50, mul1 = 1, mul2 = 10)
#' data <- afmDetachPoint(data, width = 50, mul1 = 1, mul2 = 10)
#' data <- afmBaselineCorrection(data)
#' data <- afmZeroPointSlope(data)
#' data <- afmIndentation(data)
#' data <- afmYoungModulus(data, thickness = 2e-7, params = list(alpha = 22))
#' data <- afmExpDecay(data, plt = FALSE)
#' data <- afmAdhesionEnergy(data, mul = 7)
#' 
#' # Extract the values of the parameters obtained in the analysis
#' afmExpParams <- afmExtract(data, opt.param = "type")
#' 
#' # Plotting the Young's Modulus
#' afmExpParams[[1]] %>% ggplot(aes(x = type, y = YM)) + geom_boxplot() 
#' ylab("Young's Modulus  (Pa)")
#' }
#' @export
#' 
afmExtract <- function(afmexperiment, params = list("YM", "AE", "ED"), opt.param = NULL){
  if (!is.afmexperiment(afmexperiment)){
   stop("Data should be of afmexperiment class!") 
  }
  extractedData <- data.frame(curve = names(afmexperiment))
  if (!is.null(opt.param)){
    extractedData[,eval(quote(opt.param))] <- sapply(afmexperiment, function(x) 
      get(opt.param, x$params))
  }
  if ("YM" %in% params){
    YM <- lapply(afmexperiment, function(x){ YM <- get("YoungModulus",
                                                       get("YoungModulus",x))})
    YM <- as.data.frame(do.call(rbind, YM), rownames = NULL)
    colnames(YM) <- "YM"
    extractedData <- cbind(extractedData, YM)
    row.names(extractedData) <- NULL
  }
  if ("AE" %in% params){
    AE <- lapply(afmexperiment, function(x){AE <- get("Energies", get("AdhEner",x))})
    AE <- as.data.frame(do.call(rbind, AE), rownames = NULL)
    extractedData <- cbind(extractedData, AE)
    row.names(extractedData) <- NULL
  }
  if ("ED" %in% params){
    expDecay <- lapply(afmexperiment, function(x){
      temp <- as.data.frame(coefficients(
        summary(x$ExpFit$expdecayModel))[,1:2])
      temp$parameter <- rownames(temp)
      if (!is.null(opt.param)){
        temp[,eval(quote(opt.param))]  <- get(opt.param, x$params)
      }
      return(temp)
    })
    expDecay <- do.call("rbind", expDecay)
    expDecay$curve <- as.factor(grep("force",
                                    unlist(strsplit(rownames(expDecay), ".txt.")),
                                    value = T))
    expDecay$parameter <- as.factor(expDecay$parameter)
    expDecayOrdered <- data.frame(curve = c(),Estimate = c(), StdError = c(), 
                                  parameter = c())
    if (!is.null(opt.param)){
      expDecayOrdered[,eval(quote(opt.param))] <- c()
    }
    call <- afmexperiment[[1]]$ExpFit$expdecayModel$call
    if (any(grepl("tau2", as.character(call)))){
      for (cv in levels(expDecay$curve)) {
        temp <- subset(expDecay, curve == cv)
        tau1 <- temp$Estimate[temp$parameter == "tau1"]
        tau2 <- temp$Estimate[temp$parameter == "tau2"]
        a1 <- temp$Estimate[temp$parameter == "a1"]
        a2 <- temp$Estimate[temp$parameter == "a2"]
        a0 <- temp$Estimate[temp$parameter == "a0"]
        tau1se <- temp$'Std. Error'[temp$parameter == "tau1"]
        tau2se <- temp$'Std. Error'[temp$parameter == "tau2"]
        a1se <- temp$'Std. Error'[temp$parameter == "a1"]
        a2se <- temp$'Std. Error'[temp$parameter == "a2"]
        a0se <- temp$'Std. Error'[temp$parameter == "a0"]
        if (!is.null(opt.param)) {
          op <- get(opt.param, temp)
        }
        if (tau1 > tau2) {
          temptau <- tau1
          tempc <- a1
          temptause <- tau1se
          tempcse <- a1se
          tau1 <- tau2
          tau2 <- temptau
          a1 <- a2
          a2 <- tempc
          tau1se <- tau2se
          tau2se <- temptause
          a1se <- a2se
          a2se <- tempcse
        }
        estimates <- c(a0, a1, tau1, a2, tau2)
        StdError <- c(a0se, a1se, tau1se, a2se, tau2se)
        parameter <- c("a0", "a1", "tau1", "a2", "tau2")
        if (is.null(opt.param)) {
          expDecayOrdered <- rbind(
            expDecayOrdered,
            data.frame(
              curve = cv,
              Estimate = estimates,
              StdError = StdError,
              parameter = parameter
            )
          )
        } else{
          tempdf <-
            data.frame(
              curve = cv,
              Estimate = estimates,
              StdError = StdError,
              parameter = parameter
            )
          tempdf[, eval(quote(opt.param))] <- op
          expDecayOrdered <- rbind(expDecayOrdered, tempdf)
        }
      }
    }else{
      for (cv in levels(expDecay$curve)) {
        temp <- subset(expDecay, curve == cv)
        tau1 <- temp$Estimate[temp$parameter == "tau1"]
        a1 <- temp$Estimate[temp$parameter == "a1"]
        a0 <- temp$Estimate[temp$parameter == "a0"]
        tau1se <- temp$'Std. Error'[temp$parameter == "tau1"]
        a1se <- temp$'Std. Error'[temp$parameter == "a1"]
        a0se <- temp$'Std. Error'[temp$parameter == "a0"]
        if (!is.null(opt.param)) {
          op <- get(opt.param, temp)
        }
        estimates <- c(a0, a1, tau1)
        StdError <- c(a0se, a1se, tau1se)
        parameter <- c("a0", "a1", "tau1")
        if (is.null(opt.param)) {
          expDecayOrdered <- rbind(
            expDecayOrdered,
            data.frame(
              curve = cv,
              Estimate = estimates,
              StdError = StdError,
              parameter = parameter
            )
          )
        } else{
          tempdf <-
            data.frame(
              curve = cv,
              Estimate = estimates,
              StdError = StdError,
              parameter = parameter
            )
          tempdf[, eval(quote(opt.param))] <- op
          expDecayOrdered <- rbind(expDecayOrdered, tempdf)
        }
      }
    }
    extractedData <- list(General = extractedData, ExpDecay = expDecayOrdered)
  }
  return(extractedData)
}
