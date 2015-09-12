#' Read an ascii JPK file.
#'
#' Reads an ascii JPK file with one to three headers. 
#'
#' @param filename String with the name of the jpk file
#'
#' @return A list containing a field 'data' which is a data frame
#'
#' @examples
#' afmReadJPK('testJPKfile.txt')
#'

afmReadJPK <- function(filename){
  print(filename)
  fullData <- readLines(filename)
  headerLines <- grep("#",fullData)
  numberOfHeader <- sum(diff(headerLines)!=1)+1
  N <- length(fullData)
  Nhead <- length(headerLines)
  headerStarts <- c(1,headerLines[which(diff(headerLines)!=1)+1])
  headerEnds <- c(headerLines[which(diff(headerLines)!=1)],headerLines[Nhead])
  approach <- fullData[(headerEnds[1]+1):(headerStarts[2]-2)]
  ncolumns <- length(unlist(strsplit(approach[1]," ")))
  approach <- matrix(as.numeric(unlist(strsplit(approach," "))),
                     ncol = ncolumns,
                     byrow = TRUE)
  NcolumnNames <- grep("fancy",fullData)[1]
  columnNames <- fullData[NcolumnNames]
  Fcol <- grep("Vertical" , unlist(strsplit(columnNames," \"")))-1
  Zcol <- grep("smoothed)" , unlist(strsplit(columnNames," \"")))-1
  tcol <- grep("Series", unlist(strsplit(columnNames," \"")))-1
  cnames <- c("Z","F","t")
  approach <- data.frame(Z = approach[,Zcol], 
                         Force = approach[,Fcol],
                         Time = approach[,tcol],
                         Segment = "approach" )  
  if (numberOfHeader==1){
#    retract <- NULL
#     yrange_min = min(approach$F)
#     yrange_max = max(approach$F)
#     plot(F~Z,data=approach,type="l",ylim =c(yrange_min,yrange_max) )
#     legend(x = "bottomright",legend = c("Approach"),
#            lty=1,col=c("black"))
    afmExperiment <- approach
  }
  else if(numberOfHeader==2)
  {
    retract <- fullData[(headerEnds[2]+1):N]
  
    retract <-  matrix(as.numeric(unlist(strsplit(retract," "))),
                     ncol = ncolumns,
                     byrow=TRUE)
  
    retract <- data.frame(Z = retract[,Zcol], 
                          Force = retract[,Fcol],
                          Time = retract[,tcol],
                          Segment = "retract")
#     yrange_min = min(min(approach$F),min(retract$F))
#     yrange_max = max(max(approach$F),max(retract$F))
#     plot(F~Z,data=approach,type="l",ylim =c(yrange_min,yrange_max) )
#     lines(F~Z,data=retract,col="red")
#     legend(x = "bottomright",legend = c("Approach","Retract"),
#            lty=1,col=c("black","red"))
    afmExperiment <- rbind(approach,retract)
    }
  else{
    contact <- fullData[(headerEnds[2]+1):(headerStarts[3]-2)]
    contact <-  matrix(as.numeric(unlist(strsplit(contact," "))),
                       ncol = ncolumns,
                       byrow=TRUE)
    retract <- fullData[(headerEnds[3]+1):N]
    retract <-  matrix(as.numeric(unlist(strsplit(retract," "))),
                       ncol = ncolumns,
                       byrow=TRUE)
    contact <- data.frame(Z = contact[,Zcol],
                          Force = contact[,Fcol],
                          Time = contact[,tcol],
                          Segment = "contact")
    retract <- data.frame(Z = retract[,Zcol],
                          Force = retract[,Fcol],
                          Time = retract[,tcol], 
                          Segment = "retract")
    afmExperiment <- rbind(approach, contact, retract)
  }
  afmExperiment$Segment <- factor(afmExperiment$Segment)
  return(list(data = afmExperiment))
}
