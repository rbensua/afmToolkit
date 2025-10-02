#' @title Extract computed parameters from an afmexperiment
#'   
#' @description Extracts some parameters from an afmexperiment for an easy further
#' analysis.
#' @param afmexperiment Data of afmexperiment class.
#' @param params List of parameters to extract from the data.
#' @param opt.param Optional parameter or factor in the params field of the afmdata list 
#'   to add to the data extraction.
#' @param forces When one of the parameters to extract is the "Indentation Forces" ("IF"), 
#' this numeric vector provides the forces for which obtain the corresponding force. 
#' If \code{forces} is NULL (default), then the indentation at maximum force is computed. 
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
afmExtract <- function(afmexperiment, 
                       params = list("YM", "AE", "RE","HY","AF","IF"), 
                       opt.param = NULL, 
                       forces = NULL){
  
  Segment <- NULL #For CRAN compatibility
  if (!is.afmexperiment(afmexperiment)){
   stop("Data should be of afmexperiment class!") 
  }
  extractedData <- data.frame(curve = names(afmexperiment))

## Extract optional parameters --------------
  
    if (!is.null(opt.param)){
    extractedData <- cbind(extractedData,do.call("rbind",lapply(afmexperiment, function(x) 
    as.data.frame(lapply(opt.param,function(p) get(p,x$params)), col.names = opt.param))))
    rownames(extractedData) <- NULL
  }
  
## Young modulus ---------------
  
  
  if ("YM" %in% params){
    YM <- lapply(afmexperiment, function(x){ YM <- get("YoungModulus",
                                                       get("YoungModulus",x))})
    YM <- as.data.frame(do.call(rbind, YM), rownames = NULL)
    colnames(YM) <- "YM"
    R2 <- sapply(afmexperiment, function(x){
      temp <- summary(x$YoungModulus$fitYM)
      return(temp$r.squared)
    })
    R2 <- data.frame(r.squared = R2, row.names = NULL)
    extractedData <- cbind(extractedData, YM,R2)
    row.names(extractedData) <- NULL
  }
  
## Adhesion energy --------------
  if ("AE" %in% params){
    AE <- lapply(afmexperiment, function(x){AE <- get("Energies", get("AdhEner",x))})
    AE <- as.data.frame(do.call(rbind, AE), rownames = NULL)
    extractedData <- cbind(extractedData, AE)
    row.names(extractedData) <- NULL
  }
## Relaxation parameters ---------------- 
  if ("RE" %in% params){
    if (!is.afmmulti(afmexperiment[[1]])){
      # Normal case ------------
      extractedData <- list(General = extractedData, 
                            Relax = relax(afmexperiment, opt.param = opt.param))
    } else{
      # Multi-indentation experiment
      
      numpauses <-unlist(lapply(afmexperiment, function(z) length(z$Relax$relaxModel)))
      relaxList <- list()
      for (i in 1:max(numpauses)){
        pp <- afmexperiment[which(numpauses >= i )]
        tmplist <- lapply(pp, function(X) list(params = X$params, Relax = list(model = X$Relax$model, relaxModel = X$Relax$relaxModel[[i]])))
        relaxList[[i]] <- cbind(relax(tmplist, opt.param = opt.param), data.frame(Segment = paste0("pause",i)))
      }
      extractedData <- list(General = extractedData, 
                            Relax = do.call(rbind,relaxList))
    }
  }

## Hysteresis ------------------    
  if ("HY" %in% params){
    HY <-as.data.frame(do.call(rbind,lapply(afmexperiment, 
                                            function(x) c(x$Hysteresis$Hysteresis,
                                                          x$Hysteresis$Hyst_ratio) )), 
                       rownames = NULL)
    colnames(HY) <- c("Hysteresis", "Hyst. ratio")
    extractedData <- cbind(extractedData, HY)
    row.names(extractedData) <- NULL
  }
# Adhesion Pos...? -------------------  
  if ("AF" %in% params){
    adhPos <- sapply(afmexperiment, function(x){
      ret <- subset(x$data, Segment == "retract")
      
      forceMinidx <- which.min(ret$ForceCorrected)
      adhPos <- ret$Indentation[forceMinidx] - ret$Indentation[1]
    })
    forceMin <- sapply(afmexperiment, function(x){
      ret <- subset(x$data, Segment == "retract")
      forceMin <- abs(min(ret$ForceCorrected))
    })
    dfres <- data.frame(Min.Force = forceMin, Adh.Pos = adhPos, row.names = NULL)
    extractedData <- cbind(extractedData, dfres)
  }
# Indentation Forces ----------------------  
  if ("IF" %in% params){
    
    indforces <- sapply(afmexperiment, function(x){
      app <- subset(x$data, Segment == "approach")
      if(is.null(forces)){
        forces <- max(app$ForceCorrected)
      }
      z0 <- x$Slopes$Z0Point
      idx <- sapply(seq_along(forces), function(i) max(which(app$ForceCorrected <= forces[i])))
      
      indforces <- abs(app$Indentation[idx])
      return(indforces)
    })
    A <- matrix(indforces, nrow = length(afmexperiment), ncol = length(forces), 
                byrow = TRUE,
                dimnames = list(names(afmexperiment), paste("Ind",seq_along(forces), sep = ".")))
    extractedData <- cbind(extractedData, data.frame(A, row.names = NULL))
  }
    
  
  
  return(extractedData)
}

relax <- function(afmlist, opt.param){
  nparam <- length(coefficients(afmlist[[1]]$Relax$relaxModel)) # Number of parameters of the model
  # Exponential decay models ----------------------
  if (afmlist[[1]]$Relax$model == "exp") {
    relax <- lapply(afmlist, function(x){
      temp <- as.data.frame(coefficients(
        summary(x$Relax$relaxModel))[,1:2])
      temp$parameter <- rownames(temp)
      if (!is.null(opt.param)){
        temp[,eval(quote(opt.param))]  <- get(opt.param, x$params)
      }
      return(temp)
    })
    relax <- do.call("rbind", relax)
    relax$curve <- rep(names(afmlist), each = nparam)
      # as.factor(grep(curveid,
      #                             unlist(strsplit(rownames(relax), ".txt.")),
      #                             value = T))
    relax$parameter <- as.factor(relax$parameter)
    relaxOrdered <- data.frame(curve = c(),Estimate = c(), StdError = c(), 
                               parameter = c())
    if (!is.null(opt.param)){
      relaxOrdered[,eval(quote(opt.param))] <- c()
    }
    call <- afmlist[[1]]$Relax$relaxModel$call
    type <- ifelse(any(grepl("tau", as.character(call))),"CH","CF")
    
    if(type == "CH"){
      if (any(grepl("tau2", as.character(call)))){
        for (cv in unique(relax$curve)) {
          temp <- subset(relax, curve == cv)
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
            relaxOrdered <- rbind(
              relaxOrdered,
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
            relaxOrdered <- rbind(relaxOrdered, tempdf)
          }
        }
      }else{
        for (cv in unique(relax$curve)) {
          temp <- subset(relax, curve == cv)
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
            relaxOrdered <- rbind(
              relaxOrdered,
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
            relaxOrdered <- rbind(relaxOrdered, tempdf)
          }
        }
      }
    }else{
      
      if (any(grepl("x2", as.character(call)))){
        for (cv in unique(relax$curve)) {
          temp <- subset(relax, curve == cv)
          x1 <- temp$Estimate[temp$parameter == "x1"]
          x2 <- temp$Estimate[temp$parameter == "x2"]
          c1 <- temp$Estimate[temp$parameter == "c1"]
          c2 <- temp$Estimate[temp$parameter == "c2"]
          c0 <- temp$Estimate[temp$parameter == "c0"]
          x1se <- temp$'Std. Error'[temp$parameter == "x1"]
          x2se <- temp$'Std. Error'[temp$parameter == "x2"]
          c1se <- temp$'Std. Error'[temp$parameter == "c1"]
          c2se <- temp$'Std. Error'[temp$parameter == "c2"]
          c0se <- temp$'Std. Error'[temp$parameter == "c0"]
          if (!is.null(opt.param)) {
            op <- get(opt.param, temp)
          }
          if (x1 > x2) {
            tempx <- x1
            tempc <- c1
            tempxse <- x1se
            tempcse <- c1se
            x1 <- x2
            x2 <- tempx
            c1 <- c2
            c2 <- tempc
            x1se <- x2se
            x2se <- tempxse
            c1se <- c2se
            c2se <- tempcse
          }
          estimates <- c(c0, c1, x1, c2, x2)
          StdError <- c(c0se, c1se, x1se, c2se, x2se)
          parameter <- c("c0", "c1", "x1", "c2", "x2")
          if (is.null(opt.param)) {
            relaxOrdered <- rbind(
              relaxOrdered,
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
            relaxOrdered <- rbind(relaxOrdered, tempdf)
          }
        }
      }else{
        for (cv in unique(relax$curve)) {
          temp <- subset(relax, curve == cv)
          x1 <- temp$Estimate[temp$parameter == "x1"]
          c1 <- temp$Estimate[temp$parameter == "c1"]
          c0 <- temp$Estimate[temp$parameter == "c0"]
          x1se <- temp$'Std. Error'[temp$parameter == "x1"]
          c1se <- temp$'Std. Error'[temp$parameter == "c1"]
          c0se <- temp$'Std. Error'[temp$parameter == "c0"]
          if (!is.null(opt.param)) {
            op <- get(opt.param, temp)
          }
          estimates <- c(c0, c1, x1)
          StdError <- c(c0se, c1se, x1se)
          parameter <- c("c0", "c1", "x1")
          if (is.null(opt.param)) {
            relaxOrdered <- rbind(
              relaxOrdered,
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
            relaxOrdered <- rbind(relaxOrdered, tempdf)
          }
        }
      }
      
    }
    return(relaxOrdered)
    #extractedData <- list(General = extractedData, Relax = relaxOrdered)
  } else { 
    # Power law decay models ---------------------
    relax <- lapply(afmlist, function(x){
      temp <- as.data.frame(coefficients(
        summary(x$Relax$relaxModel))[,1:2])
      temp[1,1] <- exp(temp[1,1])
      temp[1,2] <- temp[1,1]*temp[1,2]
      temp$parameter <- c("A","beta")
      if (!is.null(opt.param)){
        dfparam <- as.data.frame(lapply(opt.param, function(p) get(p,x$params)),
                                 col.names = opt.param)
        temp <- cbind(temp,dfparam)
      }
      return(temp)
    })
    relax <- do.call("rbind", relax)
    relax$curve <- rep(names(afmlist), each = nparam)
      # as.factor(grep(curveid,
      #                             unlist(strsplit(rownames(relax), ".txt.")),
      #                             value = T))
    
    rownames(relax) <- NULL
    relax <- relax[,c(ncol(relax), 1:(ncol(relax)-1))]
    #extractedData <- list(General = extractedData, Relax = relax)
    return(relax)
  }
}



