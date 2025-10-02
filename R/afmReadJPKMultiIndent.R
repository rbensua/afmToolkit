#' @title Read Nanowizard JPK ascii file containing a multi-indentation experiment
#'   
#' @description Read an ascii JPK file containing several approach - pause 
#' segments and only one retract segment.
#' 
#'
#' @param filename String with the name of the jpk file.
#' @param path Path to the folder where the file is.
#' @param FColStr String with a pattern identifying the Force column.
#' @param ZColStr String with a pattern identifying the Z column.
#' @param tColStr String with a pattern identifying the Time column.
#' @param silent Logical value. If TRUE suppresses the message indicating the name of the curve being processed
#' (useful for batch-processing large number of curves). Defaults to FALSE  
#' 
#' @return A afmdata structure list containing a field 'data' which is a data frame with 
#' variables Force, Z, Time (if aplicable) and Segment ("approach", "retract" and/or "pause") and 
#' a field 'params' which is a list with the fields 'curvename' and 'SpringConstant'.
#'   
#' @examples
#' data <- afmReadJPKMultiIndent("force-save-JPK-multiIndent.txt.gz",path = path.package("afmToolkit"))
#' # Processing
#' data <- afmContactPoint(data, width = 100, mul1 = 1, mul2 = 10)
#' data <- afmDetachPoint(data, width = 100, mul1 = 1, mul2 = 10)
#' data <- afmBaselineCorrection(data)
#' data <- afmZeroPointSlope(data, fstar = 1)
#' data <- afmIndentation(data)
#' data <- afmRelax(data, model = "power", type = "CH", tmax = 5, plt = TRUE) 
#' 
#' @importFrom stats coef lm
#' @export

afmReadJPKMultiIndent <-
  function(filename, path = "",
           FColStr = "Vertical Deflection",
           ZColStr = "Height (measured & smoothed)",
           tColStr = "Segment Time", silent = FALSE) {
    Segment <- Time <- Z <- NULL # For CRAN compatibility
    if (path == ""){
      fullfilename <- filename
    } else{
      fullfilename <- file.path(path,filename)
    }
    fullData <- readLines(fullfilename, encoding = "latin1")
    
    fullData <- fullData[which(fullData!="")]
    headerLines <- grep("#", fullData)
    
    # Obtaining the spring constant
    springLine <- grep("spring", fullData, value = T)[1]
    SpringConstant <- as.numeric(unlist(strsplit(springLine, ":"))[2])
    params = list(SpringConstant = SpringConstant, curvename = filename)
    
    
    # Obtaining the number of headers
    numberOfHeader <- as.integer(sum(diff(headerLines) != 1) + 1)
    N <- length(fullData)
    Nhead <- length(headerLines)
    
    headerStarts <- c(1, headerLines[which(diff(headerLines) != 1) + 1])
    headerEnds <-
      c(headerLines[which(diff(headerLines) != 1)], headerLines[Nhead])
    
    
    ncolumns <- length(unlist(strsplit(fullData[headerEnds[1]+1], " ")))
    NcolumnNames <- grep("fancyNames", fullData)[1]
    columnNames <- fullData[NcolumnNames]
    
    Fcol <- grep(FColStr , unlist(strsplit(columnNames, " \""))) - 1
    Fcol <- min(Fcol)
    Zcol <- which(ZColStr == gsub("\"","",unlist(strsplit(columnNames, " \""))))-1
    tcol <- grep(tColStr, unlist(strsplit(columnNames, " \""))) - 1
    cnames <- c("Z", "F", "t")
    
    afmmulti <- data.frame(Z=double(),
                     Force = double(),
                     Time = double(),
                     Segment = character(),
                     stringsAsFactors=FALSE)
    
    for (n in 1:(numberOfHeader-1)){
      curveSegment <- fullData[(headerEnds[n] + 1):(headerStarts[n+1] - 1)]
      curveSegment <- matrix(as.numeric(unlist(strsplit(curveSegment, " "))),
                         ncol = ncolumns,
                         byrow = TRUE)
      afmmulti <- rbind(afmmulti, 
                        data.frame(
                          Z = curveSegment[, Zcol],
                          Force = curveSegment[, Fcol],
                          Time = curveSegment[, tcol],
                          Segment = ifelse(n %% 2 == 1,paste0("approach",ceiling(n/2)),paste0("pause",ceiling(n/2)))
                        ))
    }
    
    retractSegment <- fullData[(headerEnds[numberOfHeader] + 1):N]
    retractSegment <- matrix(as.numeric(unlist(strsplit(retractSegment, " "))),
                            ncol = ncolumns,
                            byrow = TRUE)
    afmmulti <- rbind(afmmulti, 
                      data.frame(
                        Z = retractSegment[, Zcol],
                        Force = retractSegment[, Fcol],
                        Time = retractSegment[, tcol],
                        Segment = "retract"
                      ))
    if(!silent){
      cat(sprintf(
        "JPK file %s loaded. %d indentations found.\n",
        filename,
        numberOfHeader
      ))
    }
    afmmulti$Segment <- factor(afmmulti$Segment)
    
    ## Finding the speeds
    if ("Time" %in% colnames(afmmulti)){
      speeds <- list()
      for (segment in levels(afmmulti$Segment)){
        curve <- subset(afmmulti, Segment == segment, select = c(Time,Z))
        speed <- coefficients(lm(Z ~ Time, curve))[2]
        speeds <- append(speeds, speed)
      }
      names(speeds) <- levels(afmmulti$Segment)
    } else{
      speeds <- NULL
    }
    
    # Finding index position xpos and ypos
    indexLine <- grep("index", fullData, value = T)[1]
    Index <- as.numeric(unlist(strsplit(indexLine, ":"))[2])
    xposLine <- grep("xPosition", fullData, value = T)[1]
    xpos <- as.numeric(unlist(strsplit(xposLine, ":"))[2])
    yposLine <- grep("yPosition", fullData, value = T)[1]
    ypos <- as.numeric(unlist(strsplit(yposLine, ":"))[2])
    
    
    params <- append(params, list(speeds = speeds, Index = Index, xpos = xpos, ypos = ypos))
    res <- structure(list(data = afmmulti, params = params), class = c("afmdata","afmmulti"))
    return(res)
  }
