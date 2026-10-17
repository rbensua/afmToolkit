#' @title Read Nanowizard JPK ascii file
#'   
#' @description Read an ascii JPK file.
#' 
#' Reads an ascii JPK file with one to three headers.
#' @param filename String with the name of the jpk file.
#' @param path Path to the folder where the file is.
#' @param FColStr String with a pattern identifying the Force column.
#' @param ZColStr String with a pattern identifying the Z column.
#' @param tColStr String with a pattern identifying the Time column.
#' @param silent Logical parameter. If TRUE it suppresses the messages regarding the loading of the JPK file (defaults to FALSE).
#'   
#' @return A afmdata structure list containing a field 'data' which is a data frame with 
#' variables Force, Z, Time (if aplicable) and Segment ("approach", "retract" and/or "pause") and 
#' a field 'params' which is a list with the fields 'curvename' and 'SpringConstant'.
#'   
#' @examples
#' data <- afmReadJPK("force-save-JPK-3h.txt.gz",path = path.package("afmToolkit"))
#' str(data)
#' @importFrom stats coef lm
#' @export

afmReadJPK <-
  function(filename, 
           path = "",
           FColStr = "Vertical",
           ZColStr = "Height (measured)",
           tColStr = "Segment Time", 
           silent = FALSE) {
    
    Segment <- Time <- Z <- NULL # For CRAN compatibility
    
    # Parse the data lines of one segment into a data frame, using the
    # "fancy" column-names line of its header to locate the Force, Z and
    # Time columns. (Internal helper, defined here to avoid repeating the
    # same parsing block for the approach, pause and retract segments.)
    parseSegment <- function(lines, columnNames, segmentName) {
      ncolumns <- length(unlist(strsplit(lines[1], " ")))
      m <- matrix(as.numeric(unlist(strsplit(lines, " "))),
                  ncol = ncolumns,
                  byrow = TRUE)
      tokens <- unlist(strsplit(columnNames, " \""))
      Fcol <- grep(FColStr, tokens) - 1
      Zcol <- which(ZColStr == gsub("\"", "", tokens)) - 1
      tcol <- grep(tColStr, tokens) - 1
      data.frame(
        Z = m[, Zcol],
        Force = m[, Fcol],
        Time = m[, tcol],
        Segment = segmentName
      )
    }
    
    # JPK sometimes stores the retract segment with decreasing Z; make Z increasing.
    ensureZIncreasing <- function(segment) {
      if (coef(lm(segment$Z ~ seq_along(segment$Z)))[2] < 0) {
        segment$Z <- rev(segment$Z)
        segment$Force <- rev(segment$Force)
      }
      segment
    }
    
    if (path == ""){
      fullfilename <- filename
    } else{
      fullfilename <- file.path(path,filename)
    }
    fullData <- readLines(fullfilename)
    fullData <- fullData[sapply(fullData, nchar) > 0]
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
    # Line index of the "fancy" column-names line of each segment header
    fancyLines <- grep("fancy", fullData)
    
    # Approach segment (always the first one)
    if (numberOfHeader > 1) {
      approachLines <- fullData[(headerEnds[1] + 1):(headerStarts[2] - 2)]
    } else {
      approachLines <- fullData[(headerEnds[1] + 1):N]
    }
    approach <- parseSegment(approachLines, fullData[fancyLines[1]], "approach")
    
    if (numberOfHeader == 1) {
      afmExperiment <- approach
    }
    else if (numberOfHeader == 2){
      retractLines <- fullData[(headerEnds[2] + 1):N]
      retract <- parseSegment(retractLines, fullData[fancyLines[2]], "retract")
      retract <- ensureZIncreasing(retract)
      afmExperiment <- rbind(approach, retract)
    }
    else{
      pauseLines <- fullData[(headerEnds[2] + 1):(headerStarts[3] - 2)]
      pause <- parseSegment(pauseLines, fullData[fancyLines[2]], "pause")
      
      if (numberOfHeader > 3){
        if(!silent){
          print("More than 3 headers found, considering only the first 3!")
        }
        retractLines <- fullData[(headerEnds[3] + 1):(headerStarts[4] - 2)]
      } else{
        retractLines <- fullData[(headerEnds[3] + 1):N]
      }
      retract <- parseSegment(retractLines, fullData[fancyLines[3]], "retract")
      retract <- ensureZIncreasing(retract)
      afmExperiment <- rbind(approach, pause, retract)
    }
    if(!silent){
      cat(sprintf(
        "JPK file %s loaded. %d headers found.\n",
        filename,
        numberOfHeader
      ))
    }
    afmExperiment$Segment <- factor(afmExperiment$Segment)
    
    ## Finding the speeds
    if ("Time" %in% colnames(afmExperiment)){
      speeds <- list()
      for (segment in levels(afmExperiment$Segment)){
        curve <- subset(afmExperiment, Segment == segment, select = c(Time,Z))
        speed <- coefficients(lm(Z ~ Time, curve))[2]
        speeds <- append(speeds, speed)
      }
      names(speeds) <- levels(afmExperiment$Segment)
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
    afmExperiment <- afmdata(data = afmExperiment, params = params)
    return(afmExperiment)
  }
