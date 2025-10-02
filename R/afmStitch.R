#' @title Stitch the files of a Bruker experiment
#'
#' @description
#' Creates an afmdata structure from the four files of a Bruker experiment containing three segments:
#' approach, hold and retract. It also requires the header file to extract the spring constant.
#' @usage afmStitch(expname)
#' @param expname Veeco experiment namea character string with the common part of the 4 files.
#' @return An \code{afmdata} class variable with the three segments of the Veeco file experiment
#'
#' @examples
#' 
#' @export
afmStitch <- function(expname){
  available_files <- list.files(pattern = ".txt")
  
  
  # final part of the curve files names
  finalpart <- c("approach.txt","hold.txt","retract.txt")
  filenames <- paste(expname,paste("seg",1:3,sep = ""),finalpart,sep  = "_")
  # header  
  header <- paste(expname, "header.txt", sep = "_")
  
  # checking if all the necessary files are there
  if (!all(c(filenames,header)%in% available_files)){
    error_message <- paste("Some file parts are missing for experiment",expname, sep = " ")
    stop(error_message)
  }
  
  springLine <- grep("spring",readLines(header), value = TRUE, ignore.case = TRUE)[1]
  SpringConstant <-
    as.numeric(sub("\"", "", unlist(strsplit(springLine, ":"))[2]))

  
  #column names  
  cnames <- c("Time","Force","Z","Force2","Z2")
  
  # set the afmdata parameter list
  params <- list(SpringConstant = SpringConstant, curvename = expname)
  
  # Reading the three segments
  seg1 <- read.table(filenames[1], sep = "\t", skip = 11)
  seg2 <- read.table(filenames[2], sep = "\t", skip = 11)
  seg3 <- read.table(filenames[3], sep = "\t", skip = 11)
  
  # setting the column names
  colnames(seg1) <- cnames
  colnames(seg2) <- cnames
  colnames(seg3) <- cnames
  
  # getting rid of the last two columns (which seem to be mere repetions)
  seg1 <- seg1[ ,1:3]
  seg2 <- seg2[ ,1:3]
  seg3 <- seg3[ ,1:3]
  
  # Creating the Segment identification column
  seg1$Segment = "approach"
  seg2$Segment = "pause"
  seg3$Segment = "retract"
  
  # Binding the three segments into a single data.frame 
  # (which will be the "data" field of the afmdata structure)
  dataFull <- rbind(seg1,seg2,seg3)
  dataFull$Z <- dataFull$Z*1e-6
  dataFull$Segment <- as.factor(dataFull$Segment)
  
  # Creating the afmdata structure
  dataFull <- afmdata(dataFull, params = params)
  
  return(dataFull)
}