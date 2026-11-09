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
#' @details The Force, Z and Time columns of each segment are located in the
#' "fancy names" line of the JPK header: \code{FColStr}, \code{ZColStr} and
#' \code{tColStr} are first tried as exact matches, then as substrings and
#' finally ignoring spaces and parentheses, so minor variations of the channel
#' names in the export (e.g. "Height (measured)" matching "Height (measured &
#' smoothed)") are still recognized. When several channels match the force
#' pattern, deflection channels are preferred. If no unique column can be
#' identified, an informative error is thrown. Files with a single segment
#' take their segment name ("approach", "retract" or "pause") from the header,
#' defaulting to "approach" when the header does not declare one.
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
    
    # Locate a column of a data block among the "fancy" channel names: an
    # exact match is preferred over a substring match, which is preferred
    # over a match ignoring spaces and parentheses; when several channels
    # match the pattern, prefer (if not NULL) tries to disambiguate. Returns
    # the index of the column in the data matrix, or fails with an
    # informative error.
    findCol <- function(pattern, tokens, prefer = NULL) {
      hits <- which(tokens == pattern)
      if (length(hits) == 0) {
        hits <- grep(pattern, tokens, fixed = TRUE)
      }
      if (length(hits) == 0) {
        # Last resort: compare ignoring spaces and parentheses, so e.g.
        # "Height (measured)" also matches "Height (measured & smoothed)".
        hits <- grep(gsub("[ ()]", "", pattern),
                     gsub("[ ()]", "", tokens),
                     fixed = TRUE)
      }
      if (length(hits) > 1 && !is.null(prefer)) {
        preferred <- grep(prefer, tokens[hits], fixed = TRUE)
        if (length(preferred) == 1) {
          hits <- hits[preferred]
        }
      }
      if (length(hits) != 1) {
        stop("afmReadJPK: could not uniquely identify the column matching '",
             pattern, "'. Available columns: ",
             paste(tokens, collapse = ", "),
             ". Please adjust FColStr, ZColStr or tColStr accordingly.",
             call. = FALSE)
      }
      hits
    }
    
    # Parse the data lines of one segment into a data frame, using the
    # "fancy" column-names line of its header to locate the Force, Z and
    # Time columns. (Internal helper, defined here to avoid repeating the
    # same parsing block for the approach, pause and retract segments.)
    parseSegment <- function(lines, columnNames, segmentName) {
      ncolumns <- length(unlist(strsplit(lines[1], " ")))
      m <- matrix(as.numeric(unlist(strsplit(lines, " "))),
                  ncol = ncolumns,
                  byrow = TRUE)
      tokens <- gsub("\"", "", unlist(strsplit(columnNames, " \"")))
      tokens <- tokens[!grepl("^#", tokens)]
      Fcol <- findCol(FColStr, tokens, prefer = "Deflection")
      Zcol <- findCol(ZColStr, tokens)
      tcol <- findCol(tColStr, tokens)
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
    
    # Segment type declared in the header lines ("approach", "retract" or
    # "pause"); NA when the header does not declare one.
    segmentInHeader <- function(lines) {
      segLine <- grep("# segment:", lines, fixed = TRUE, value = TRUE)
      if (length(segLine) == 0) {
        return(NA_character_)
      }
      trimws(sub("# segment:", "", segLine[1]))
    }
    
    if (path == ""){
      fullfilename <- filename
    } else{
      fullfilename <- file.path(path,filename)
    }
    fullData <- readLines(fullfilename)
    fullData <- fullData[sapply(fullData, nchar) > 0]
    headerLines <- grep("#", fullData)
    
    # Obtaining the spring constant (absent in some exports, e.g. map files)
    springLine <- grep("spring", fullData, value = TRUE)
    if (length(springLine) == 0) {
      SpringConstant <- NA_real_
      if (!silent) {
        warning("No spring constant found in the JPK file.")
      }
    } else {
      SpringConstant <- as.numeric(unlist(strsplit(springLine[1], ":"))[2])
    }
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
    # Single-segment files are not always approach curves: take the segment
    # name declared in the header when present.
    firstSegmentName <- "approach"
    if (numberOfHeader == 1) {
      declaredSegment <- segmentInHeader(fullData[1:headerEnds[1]])
      if (!is.na(declaredSegment)) {
        firstSegmentName <- declaredSegment
      }
    }
    approach <- parseSegment(approachLines, fullData[fancyLines[1]], firstSegmentName)
    
    if (numberOfHeader == 1) {
      afmExperiment <- approach
      if (firstSegmentName == "retract") {
        afmExperiment <- ensureZIncreasing(afmExperiment)
      }
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
    
    # Finding index position xpos and ypos (these keys are absent from some
    # exports, e.g. single-curve map files; the values then stay NA)
    numFromHeader <- function(pattern) {
      line <- grep(pattern, fullData, value = TRUE)
      if (length(line) == 0) {
        return(NA_real_)
      }
      as.numeric(unlist(strsplit(line[1], ":"))[2])
    }
    Index <- numFromHeader("index")
    xpos <- numFromHeader("xPosition")
    ypos <- numFromHeader("yPosition")
    
    
    params <- append(params, list(speeds = speeds, Index = Index, xpos = xpos, ypos = ypos))
    afmExperiment <- afmdata(data = afmExperiment, params = params)
    return(afmExperiment)
  }
