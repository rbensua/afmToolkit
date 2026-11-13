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
#' "retract" (last block). When several header blocks declare the same
#' segment name (e.g. several "modulation" segments), the data frame also
#' carries a \code{SegmentID} column labelling every occurrence uniquely
#' in file order ("modulation_1", "modulation_2", ...); segment names
#' appearing only once keep their plain name, and files without repeated
#' segment names have no \code{SegmentID} column. The data of each segment
#' is located through its
#' "# segmentIndex:" key when the file declares one per segment header (JPK
#' 5.0 and 6.1 exports do, possibly writing all the header blocks
#' consecutively, without blank lines in between), and through the runs of
#' "#" lines separated by blank lines otherwise; only the numeric rows of
#' each segment are kept, so the "#" separator and key lines are dropped
#' whatever the export layout. The Force, Z and Time columns of each segment
#' are located in the "fancy names" line of its header: \code{FColStr},
#' \code{ZColStr} and \code{tColStr} are first tried as exact matches, then as
#' substrings and finally ignoring spaces and parentheses, so minor variations
#' of the channel names in the export (e.g. "Height (measured)" also matching
#' "Height (measured & smoothed)") are still recognized. When several channels
#' match the force pattern, deflection channels are preferred. If no unique
#' column can be identified, an informative error is thrown.
#' 
#' @return A afmdata structure list containing a field 'data' which is a data frame with 
#' variables Force, Z, Time (if aplicable) and Segment (one level per segment of 
#' the file, see Details), a SegmentID column when the file has several
#' segments of the same type (see Details) and a field 'params' which is a list with the fields 
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
      if (length(lines) == 0) {
        return(data.frame(Z = numeric(0), Force = numeric(0),
                          Time = numeric(0), Segment = segmentName))
      }
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
    N <- length(fullData)
    # Line index of the "fancyNames" column-names line of each segment
    # header (the plain "fancy" pattern would also match the
    # "prefix-fancy"/"suffix-fancy" keys written by JPK 5.0 exports)
    fancyLines <- grep("fancyNames", fullData, fixed = TRUE)
    if (length(fancyLines) == 0) {
      stop("afmReadJPK: no segment headers (fancyNames key) found in the file.",
           call. = FALSE)
    }
    # Line index of the "# segmentIndex:" key of each segment header
    segIdxLines <- grep("# segmentIndex:", fullData, fixed = TRUE)
    
    # Data range of each segment. When the file declares one
    # "# segmentIndex:" key per segment header (JPK 5.0 and 6.1 exports do,
    # possibly writing all the header blocks consecutively, without blank
    # lines in between), the data of a segment runs from the line after its
    # fancyNames line to the line before the next "# segmentIndex:" key;
    # otherwise the runs of "#" lines separated by blank lines delimit the
    # header blocks, as in the original reader.
    if (length(segIdxLines) == length(fancyLines)) {
      numberOfHeader <- length(segIdxLines)
      dataStarts <- fancyLines + 1
      dataEnds <- c(segIdxLines[-1] - 1, N)
      nameBlockStarts <- segIdxLines
    } else {
      numberOfHeader <- as.integer(sum(diff(headerLines) != 1) + 1)
      Nhead <- length(headerLines)
      headerStarts <- c(1, headerLines[which(diff(headerLines) != 1) + 1])
      headerEnds <-
        c(headerLines[which(diff(headerLines) != 1)], headerLines[Nhead])
      dataStarts <- headerEnds + 1
      dataEnds <- c(headerStarts[-1] - 2, N)
      nameBlockStarts <- headerStarts
    }
    if (numberOfHeader != length(fancyLines)) {
      stop("afmReadJPK: found ", numberOfHeader, " header blocks for ",
           length(fancyLines), " segment headers; unsupported file layout.",
           call. = FALSE)
    }
    
    # Segment name of each header block: taken from the "# segment:" key of
    # the block ("extend" is renamed to "approach", any other declared name
    # is kept as is); blocks without the key fall back to the positional
    # names approach/pause/retract.
    segmentNames <- vapply(seq_len(numberOfHeader), function(i) {
      declared <- segmentInHeader(fullData[nameBlockStarts[i]:fancyLines[i]])
      if (!is.na(declared)) {
        declared
      } else if (i == 1) {
        "approach"
      } else if (i == numberOfHeader) {
        "retract"
      } else {
        "pause"
      }
    }, character(1))
    # Occurrence index of each header block within its segment name: when
    # several blocks declare the same name (e.g. several "modulation"
    # segments), each one gets a unique SegmentID "name_1", "name_2", ...
    # in file order; names appearing only once keep their plain name.
    nameCounts <- table(segmentNames)[segmentNames]
    occurrence <- ave(seq_len(numberOfHeader), segmentNames, FUN = seq_along)
    segmentIDs <- ifelse(nameCounts > 1,
                         paste0(segmentNames, "_", occurrence),
                         segmentNames)
    
    # Parse the data lines of every segment: only the numeric rows are
    # kept, so the "#"-key and separator lines inside the range are dropped
    # whatever the export layout. Retract segments get the Z-increasing
    # normalisation.
    segments <- lapply(seq_len(numberOfHeader), function(i) {
      if (dataEnds[i] >= dataStarts[i]) {
        segmentLines <- fullData[dataStarts[i]:dataEnds[i]]
      } else {
        segmentLines <- character(0)
      }
      rowTokens <- strsplit(trimws(segmentLines), "\\s+")
      isData <- vapply(rowTokens, function(tk) {
        length(tk) > 0 && all(!is.na(suppressWarnings(as.numeric(tk))))
      }, logical(1))
      segment <- parseSegment(segmentLines[isData], fullData[fancyLines[i]], segmentNames[i])
      if (segmentNames[i] == "retract") {
        segment <- ensureZIncreasing(segment)
      }
      segment
    })
    afmExperiment <- do.call(rbind, segments)
    # When the file has several segments of the same type, add a SegmentID
    # column labelling every occurrence uniquely ("modulation_1",
    # "modulation_2", ...); files without repeated segment names keep the
    # plain four-column data frame.
    if (any(nameCounts > 1)) {
      afmExperiment$SegmentID <- factor(rep(segmentIDs,
                                           vapply(segments, nrow, integer(1))))
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
