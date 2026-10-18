# Internal module of summary.afmdata(): builds and prints the numeric
# (console) summary of an afmdata object, based on which analysis fields
# are present. Returns (invisibly) a list with the summary data frame
# (output) and the type of the relaxation fit ("Constant Force" /
# "Constant Height", NULL when absent), which the pause-segment plot
# module needs.
afmSummaryTable <- function(afmdata) {
  typeoffit <- NULL
  fields <- names(afmdata)
  Nsegments <- nlevels(afmdata$data$Segment)
  spring <- afmdata$params$SpringConstant
  output <- data.frame(Nsegments, spring, row.names = NULL)
  names(output) <- c("# of segments", "Spring Constant")
  if ("CP" %in% fields) {
    cp <- afmdata$CP$CP
    names(cp) <- NULL
    newnames <- c(names(output), "Contact Point")
    output <- cbind(output, cp)
    names(output) <- newnames
  }
  if ("YoungModulus" %in% names(afmdata)) {
    ym <- afmdata$YoungModulus$YoungModulus
    newnames <- c(names(output), "Young's Modulus")
    output <- cbind(output, ym)
    names(output) <- newnames
  }
  if ("Slopes" %in% names(afmdata)) {
    z0 <- afmdata$Slopes$Z0Point
    names(z0) <- NULL
    newnames <- c(names(output), "Zero force pt.")
    output <- cbind(output, z0)
    names(output) <- newnames
  }
  if ("Relax" %in% names(afmdata)) {
    if (grepl("Z ~", afmdata$Relax$relaxModel$call[2])) {
      typeoffit <- "Constant Force"
    } else {
      typeoffit <- "Constant Height"
    }
    newnames <- c(names(output), "Type of Experiment")
    output <- cbind(output, typeoffit)
    names(output) <- newnames
    Relax <- summary(afmdata$Relax$relaxModel)
    print(Relax$coefficients)
  }
  print(output)
  invisible(list(output = output, typeoffit = typeoffit))
}
