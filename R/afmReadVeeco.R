#' @title Read Bruke Nanoscope Veeco ascii file
#'   
#' @description Read an ascii Veeco file.
#' 
#' Reads an ascii Veeco file with one or two segments.
#' @param filename String with the name of the jpk file.
#' @param path Path to the folder where the file is.
#' @param FColStr String pattern identifying the Force columns (defaults to "pN")
#' @param ZColStr String pattern identifying the Z columns (defaults to "Ramp")
#' @param tColStr String pattern identifying the Time columns (defaults to "Time")
#' @param TimeCol Logical value. If TRUE (default) there is a Time column.
#' @param silent Logical parameter. If TRUE it suppresses the messages regarding the loading of the JPK file (defaults to FALSE).
#' 
#' @return A afmdata structure list containing a field 'data' which is a data frame with 
#' variables Force, Z, Time (if aplicable) and Segment ("approach", "retract" and/or "pause") and 
#' a field 'params' which is a list with the fields 'curvename' and 'SpringConstant'.
#'  
#' @examples
#' data <- afmReadVeeco("veeco_file.txt.gz",path = path.package("afmToolkit"))
#' str(data)
#' @importFrom stats coef lm
#' @export

afmReadVeeco <-
  function(filename, 
           path = ".", 
           FColStr = "pN",
           ZColStr = "Ramp",
           tColStr = "Time", 
           TimeCol = TRUE, 
           silent = FALSE) {
    if (path == "") {
      fullfilename <- filename
    } else{
      fullfilename <- file.path(path, filename)
    }
    fullData <- readLines(fullfilename)
    fullData <- fullData[sapply(fullData, nchar) > 0]
    # Obtaining the spring constant
    springLine <- grep("Spring", fullData, value = T)[1]
    SpringConstant <-
      as.numeric(sub("\"", "", unlist(strsplit(springLine, ":"))[2]))
    params = list(SpringConstant = SpringConstant, curvename = filename)
    # Obtaining the number of headers
    Data <- grep("\"", fullData, invert = TRUE, value = TRUE)
    #DataDF <- read.table(text = Data, sep ="\t", header = T)
    DataDF <- read.table(text = Data, header = T)
    cnames <- colnames(DataDF)
    # Getting the approach segment columns (Extent)
    appCols  <- grep("Ex", cnames, value = "T")
    if (TimeCol) {
      appTime <- grep(tColStr, appCols, value = "T")
      appTime <- grep(appTime, cnames)
    }
    appZ <- grep(ZColStr, appCols, value = "T")
    appForce <- grep(FColStr, appCols, value = "T")
    appZ <- grep(appZ, cnames)
    appForce <- grep(appForce, cnames)
    if (TimeCol) {
      appSegment <- data.frame(
        Time = DataDF[, appTime],
        Z = DataDF[, appZ] * 1e-9,
        Force = DataDF[, appForce] * 1e-12,
        Segment = "approach"
      )
    } else{
      appSegment <- data.frame(Z = DataDF[, appZ] * 1e-9,
                               Force = DataDF[, appForce] * 1e-12,
                               Segment = "approach")
    }
    appSegment <- appSegment[complete.cases(appSegment), ]
    
    
    # Getting the retract segment columns (Retract)
    retCols  <- grep("Rt", cnames, value = "T")
    if (TimeCol) {
      retTime <- grep(tColStr, retCols, value = "T")
      retTime <- grep(retTime, cnames)
    }
    retZ <- grep(ZColStr, retCols, value = "T")
    retForce <- grep(FColStr, retCols, value = "T")
    retZ <- grep(retZ, cnames)
    retForce <- grep(retForce, cnames)
    if(TimeCol) {
      retSegment <- data.frame(
        Time = DataDF[, retTime],
        Z = DataDF[, retZ] * 1e-9,
        Force = DataDF[, retForce] * 1e-12,
        Segment = "retract"
      )
    } else{
      retSegment <- data.frame(Z = DataDF[, retZ] * 1e-9,
                               Force = DataDF[, retForce] * 1e-12,
                               Segment = "retract")
    }
    retSegment <- retSegment[complete.cases(retSegment), ]
    napp <- nrow(appSegment)
    applinefit <- lm(Z ~ n, data = data.frame(Z = appSegment$Z, n = 1:napp))
    if (coefficients(applinefit)["n"] > 0)
    {
      retSegment$Z <- rev(retSegment$Z)
      appSegment$Z <- rev(appSegment$Z)
    }
    afmExperiment <- rbind(appSegment, retSegment)
    if (!silent) {
      cat(sprintf("Veeco file %s loaded.\n",
                  filename))
    }
    afmExperiment$Segment <- as.factor(afmExperiment$Segment)
    afmExperiment <- afmdata(afmExperiment, params = params)
    return(afmExperiment)
  }
