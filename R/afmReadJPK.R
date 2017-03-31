#' @title Read Nanowizard JPK ascii file
#'   
#' @description Read an ascii JPK file.
#' 
#' Reads an ascii JPK file with one to three headers.
#' @usage afmReadJPK(filename, path = "", FColStr = "Vertical", 
#' ZColStr = "Height (measured & smoothed)", tColStr = "Segment Time")
#' @param filename String with the name of the jpk file.
#' @param path Path to the folder where the file is.
#' @param FColStr String with a pattern identifying the Force column.
#' @param ZColStr String with a pattern identifying the Z column.
#' @param tColStr String with a pattern identifying the Time column.
#'   
#' @return A list containing a field 'data' which is a data frame
#'   
#' @examples
#' data <- afmReadJPK("force-save-JPK-3h.txt.gz",path = path.package("afmToolkit"))
#' str(data)
#' @importFrom stats coef lm
#' @export

afmReadJPK <-
  function(filename, path = "",
           FColStr = "Vertical",
           ZColStr = "Height (measured & smoothed)",
           tColStr = "Segment Time") {
    fullfilename <- file.path(path,filename)
    fullData <- readLines(fullfilename)
    fullData <- fullData[sapply(fullData, nchar) > 0]
    headerLines <- grep("#", fullData)
    
    # Obtaining the spring constant
    springLine <- grep("spring", fullData, value = T)[1]
    SpringConstant <- as.numeric(unlist(strsplit(springLine, ":"))[2])
    params = list(SpringConstant = SpringConstant, curvename = filename)
    # Obtaining the number of headers
    numberOfHeader <- as.integer(sum(diff(headerLines) != 1) + 1)
    if (numberOfHeader>3){
      stop("Currently only up to three segments are supported!")
    }
    N <- length(fullData)
    Nhead <- length(headerLines)
    
    headerStarts <- c(1, headerLines[which(diff(headerLines) != 1) + 1])
    headerEnds <-
      c(headerLines[which(diff(headerLines) != 1)], headerLines[Nhead])
    if (numberOfHeader > 1) {
      approach <- fullData[(headerEnds[1] + 1):(headerStarts[2] - 2)]
    } else if (numberOfHeader == 1) {
      approach <- fullData[(headerEnds[1] + 1):N]
    }
    ncolumns <- length(unlist(strsplit(approach[1], " ")))
    approach <- matrix(as.numeric(unlist(strsplit(approach, " "))),
                       ncol = ncolumns,
                       byrow = TRUE)
    NcolumnNames <- grep("fancy", fullData)[1]
    columnNames <- fullData[NcolumnNames]
    Fcol <- grep(FColStr , unlist(strsplit(columnNames, " \""))) - 1
#   Zcol <- grep(ZColStr , unlist(strsplit(columnNames, " \""))) - 1
    Zcol <- which(ZColStr == gsub("\"","",unlist(strsplit(columnNames, " \""))))-1
    tcol <- grep(tColStr, unlist(strsplit(columnNames, " \""))) - 1
    cnames <- c("Z", "F", "t")
    approach <- data.frame(
      Z = approach[, Zcol],
      Force = approach[, Fcol],
      Time = approach[, tcol],
      Segment = "approach"
    )
    if (numberOfHeader == 1) {
      afmExperiment <- approach
    }
    else if (numberOfHeader == 2)
    {
      retract <- fullData[(headerEnds[2] + 1):N]
      
      retract <-  matrix(as.numeric(unlist(strsplit(retract, " "))),
                         ncol = ncolumns,
                         byrow = TRUE)
      
      retract <- data.frame(
        Z = retract[, Zcol],
        Force = retract[, Fcol],
        Time = retract[, tcol],
        Segment = "retract"
      )
      if (coef(lm(retract$Z ~ seq_along(retract$Z)))[2] < 0) {
        retract$Z <- rev(retract$Z)
        retract$Force <- rev(retract$Force)
      }
      afmExperiment <- rbind(approach, retract)
    }
    else{
      pause <- fullData[(headerEnds[2] + 1):(headerStarts[3] - 2)]
      pause <-  matrix(as.numeric(unlist(strsplit(pause, " "))),
                         ncol = ncolumns,
                         byrow = TRUE)
      retract <- fullData[(headerEnds[3] + 1):N]
      retract <-  matrix(as.numeric(unlist(strsplit(retract, " "))),
                         ncol = ncolumns,
                         byrow = TRUE)
      pause <- data.frame(
        Z = pause[, Zcol],
        Force = pause[, Fcol],
        Time = pause[, tcol],
        Segment = "pause"
      )
      retract <- data.frame(
        Z = retract[, Zcol],
        Force = retract[, Fcol],
        Time = retract[, tcol],
        Segment = "retract"
      )
      if (coef(lm(retract$Z ~ seq_along(retract$Z)))[2] < 0) {
        retract$Z <- rev(retract$Z)
        retract$Force <- rev(retract$Force)
      }
      afmExperiment <- rbind(approach, pause, retract)
    }
    cat(sprintf(
      "JPK file %s loaded. %d headers found.\n",
      filename,
      numberOfHeader
    ))
    
    afmExperiment$Segment <- factor(afmExperiment$Segment)
    
    afmExperiment <- afmdata(data = afmExperiment, params = params)
    return(afmExperiment)
  }
