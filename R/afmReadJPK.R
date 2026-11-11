#' @title Read Nanowizard JPK ascii file
#'   
#' @description Read an ascii JPK file.
#' 
#' Reads an ascii JPK file with any number of segment headers.
#' @param filename String with the name of the jpk file.
#' @param path Path to the folder where the file is.
#' @param FColStr String with a pattern identifying the Force column.
#' @param ZColStr String with a pattern identifying the Z column.
#' @param tColStr String with a pattern identifying the Time column.
#' @param silent Logical parameter. If TRUE it suppresses the messages regarding the loading of the JPK file (defaults to FALSE).
#'   
#' @details All the segment headers of the file are read, and the name of each
#' segment is taken from its "# segment:" key: "extend" is renamed to
#' "approach" and any other declared name ("pause", "retract", "modulation",
#' ...) is kept as is. Header blocks without a "# segment:" key fall back to
#' the positional names "approach" (first block), "pause" (middle blocks) and
#' "retract" (last block). The Force, Z and Time columns of each segment are
#' located in the "fancy names" line of its header: \code{FColStr},
#' \code{ZColStr} and \code{tColStr} are first tried as exact matches, then as
#' substrings and finally ignoring spaces and parentheses, so minor variations
#' of the channel names in the export (e.g. "Height (measured)" also matching
#' "Height (measured & smoothed)") are still recognized. When several channels
#' match the force pattern, deflection channels are preferred. If no unique
#' column can be identified, an informative error is thrown.
#' 
#' @return A afmdata structure list containing a field 'data' which is a data frame with 
#' variables Force, Z, Time (if aplicable) and Segment (one level per segment of 
#' the file, see Details) and a field 'params' which is a list with the fields 
#' 'curvename' and 'SpringConstant'.
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
    
    # Segment type declared in the header lines; "extend" is renamed to
    # "approach" and any other declared name is kept as is. NA when the
    # header does not declare one.
    segmentInHeader <- function(lines) {
      segLine <- grep("# segment:", lines, fixed = TRUE, value = TRUE)
      if (length(segLine) == 0) {
        return(NA_character_)
      }
      segmentName <- trimws(sub("# segment:", "", segLine[1]))
      if (segmentName == "extend") {
        segmentName <- "approach"
      }
      segmentName
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
    
    # Segment name of each header block: taken from the "# segment:" key of
    # the block ("extend" is renamed to "approach", any other declared name
    # is kept as is); blocks without the key fall back to the positional
    # names approach/pause/retract.
    segmentNames <- character(numberOfHeader)
    for (i in seq_len(numberOfHeader)) {
      declared <- segmentInHeader(fullData[headerStarts[i]:headerEnds[i]])
      if (!is.na(declared)) {
        segmentNames[i] <- declared
      } else if (i == 1) {
        segmentNames[i] <- "approach"
      } else if (i == numberOfHeader) {
        segmentNames[i] <- "retract"
      } else {
        segmentNames[i] <- "pause"
      }
    }
    
    # Parse the data lines of every segment: non-final segments run up to
    # the line before the next header block and the last one runs to the
    # end of the file. Retract segments get the Z-increasing normalisation.
    segments <- lapply(seq_len(numberOfHeader), function(i) {
      if (i < numberOfHeader) {
        segmentLines <- fullData[(headerEnds[i] + 1):(headerStarts[i + 1] - 2)]
      } else {
        segmentLines <- fullData[(headerEnds[i] + 1):N]
      }
      segment <- parseSegment(segmentLines, fullData[fancyLines[i]], segmentNames[i])
      if (segmentNames[i] == "retract") {
        segment <- ensureZIncreasing(segment)
      }
      segment
    })
    afmExperiment <- do.call(rbind, segments)
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
