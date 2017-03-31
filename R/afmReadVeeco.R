#' @title Read Bruke Nanoscope Veeco ascii file
#'   
#' @description Read an ascii Veeco file.
#' 
#' Reads an ascii Veeco file with one or two segments.
#' @usage afmReadVeeco(filename, path = "")
#' @param filename String with the name of the jpk file.
#' @param path Path to the folder where the file is.
#' @return A list containing a field 'data' which is a data frame
#'   
#' @examples
#' data <- afmReadVeeco("veeco_file.txt.gz",path = path.package("afmToolkit"))
#' str(data)
#' @importFrom stats coef lm
#' @export

afmReadVeeco <-
  function(filename, path = "") {
    fullfilename <- file.path(path,filename)
    fullData <- readLines(fullfilename)
    fullData <- fullData[sapply(fullData, nchar) > 0]
    # Obtaining the spring constant
    springLine <- grep("Spring", fullData, value = T)[1]
    SpringConstant <- as.numeric(sub("\"","",unlist(strsplit(springLine, ":"))[2]))
    params = list(SpringConstant = SpringConstant, curvename = filename)
    # Obtaining the number of headers
    Data <- grep("\"",fullData,invert = TRUE, value = TRUE)
    DataDF <- read.table(text = Data, sep ="\t", header = T)
    cnames <- colnames(DataDF)
    appCols  <- grep("Ex",cnames, value = "T")
    appTime <- grep("Time",appCols, value = "T")
    appZ <- grep("Ramp", appCols, value = "T")
    appForce <- grep("pN", appCols, value = "T")
    appTime <- grep(appTime, cnames)
    appZ <- grep(appZ, cnames)
    appForce <- grep(appForce, cnames)
    appSegment <- data.frame(Time = DataDF[,appTime], 
                             Z = DataDF[,appZ]*1e-9, 
                             Force = DataDF[,appForce]*1e-12,
                             Segment = "approach")
    appSegment <- appSegment[complete.cases(appSegment),]
    
    retCols  <- grep("Rt",cnames, value = "T")
  retTime <- grep("Time",retCols, value = "T")
    retZ <- grep("Ramp", retCols, value = "T")
    retForce <- grep("pN", retCols, value = "T")
  retTime <- grep(retTime, cnames)
    retZ <- grep(retZ, cnames)
    retForce <- grep(retForce, cnames)
    retSegment <- data.frame(Time = DataDF[,retTime], 
                             Z = DataDF[,retZ]*1e-9, 
                             Force = DataDF[,retForce]*1e-12,
                             Segment = "retract")
    retSegment <- retSegment[complete.cases(retSegment),]
    if (diff(appSegment$Z)[1]>0)
    {
      retSegment$Z <- rev(retSegment$Z)
      appSegment$Z <- rev(appSegment$Z)
    }
    afmExperiment <- rbind(appSegment,retSegment)
    cat(sprintf(
      "Veeco file %s loaded.\n",
      filename
    ))
    afmExperiment$Segment <- as.factor(afmExperiment$Segment)
    afmExperiment <- afmdata(afmExperiment, params = params)
    return(afmExperiment)
  }
