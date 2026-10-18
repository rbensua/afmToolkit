# Internal module of summary.afmdata(): builds and prints the numeric
# (console) summary of an afmdata object, based on which analysis fields
# are present. Returns (invisibly) a list with the summary data frame
# (output) and the type of the relaxation fit ("Constant Force" /
# "Constant Height", NULL when absent), which the pause-segment plot
# module needs.
afmSummaryTable <- function(afmdata) {
  typeoffit <- NULL
  relax_formula <- NULL
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
    relax_formula <-
      paste(deparse(afmdata$Relax$relaxModel$call[[2]]), collapse = " ")
    if (grepl("Z ~", relax_formula)) {
      typeoffit <- "Constant Force"
    } else {
      typeoffit <- "Constant Height"
    }
    newnames <- c(names(output), "Type of Experiment")
    output <- cbind(output, typeoffit)
    names(output) <- newnames
  }
  # The data frame always has a single row, so the implicit row name
  # ("1") carries no information: print without it.
  print(output, row.names = FALSE)
  # Relaxation fit: shown as its own clearly separated section, with the
  # fitted formula (taken from the model call) and its coefficients table.
  if (!is.null(typeoffit)) {
    cat("\n", strrep("-", 60), "\n", sep = "")
    cat("Relaxation fit (model: ", afmdata$Relax$model, ", ", typeoffit, "):\n", sep = "")
    cat("  ", relax_formula, "\n\n", sep = "")
    print(summary(afmdata$Relax$relaxModel)$coefficients)
    cat(strrep("-", 60), "\n", sep = "")
  }
  invisible(list(output = output, typeoffit = typeoffit))
}
