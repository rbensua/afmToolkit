#' @title Plot an afmdata object
#' 
#' 
#' @description 
#' Plots an afmdata object.
#' 
#' @param \code{afmdata}: An object of \code{afmdata} class.
#' @param \code{vs}: The variable for the x-axis. May take the values "Time" or 
#'   "Z". It defaults to "Z", plotting thus a Force-Distance curve. If \code{vs}
#'   is set to "Time", then it plots  a Force-Time curve.
#' @param \code{segment}: The segment of the curve to be plotted. If
#'   \code{segment = "all"} then all segments of the curve are plotted. Possible
#'   values are: \code{"approach"}, \code{"contact"}, \code{"retract"} and 
#'   \code{"all"}.
#'   
#' 
#' @export
#' 

summary.afmdata <- function(afmdata,plt = TRUE){
  if (!is.afmdata(afmdata)){
    stop("Input should be an afmdata structure!")
  }
  nfields <- length(afmdata)
  fields <- names(afmdata)
  Nsegments <- nlevels(afmdata$data$Segment)
  spring <- afmdata$params$SpringConstant
  output <- data.frame(Nsegments, spring, row.names = NULL)
  names(output) <- c("# of segments", "Spring Constant")
  if ("CP" %in% fields){
    cp <- afmdata$CP$CP
    names(cp) <- NULL
    newnames <- c(names(output),"Contact Point")
    output <- cbind(output,cp)
    names(output) <- newnames
  }
  if ("YoungModulus" %in% names(afmdata)){
    ym <- afmdata$YoungModulus$YoungModulus
    newnames <- c(names(output),"Young's Modulus")
    output <- cbind(output,ym)
    names(output) <- newnames
  }
  if ("Slopes" %in% names(afmdata)){
    z0 <- afmdata$Slopes$Z0Point
    names(z0) <- NULL
    newnames <- c(names(output),"Zero force pt.")
    output <- cbind(output,z0)
    names(output) <- newnames
  }
  if ("ExpFit" %in% names(afmdata)){
    if(grepl("Force ~",afmdata$ExpFit$expdecayModel$call[2])){
      typeoffit <- "Constant Height"
    }else{
      typeoffit <- "Constant Force"
    }
    newnames <- c(names(output),"Type of Experiment")
    output <- cbind(output,typeoffit)
    names(output) <- newnames
    expfit <- summary(afmdata$ExpFit$expdecayModel)
  }
  if (plt){
    
    fancy_scientific <- function(l) {
      # turn in to character string in scientific notation
      l <- format(l, scientific = TRUE)
      l[as.numeric(l) == 0] <- "0"
      # quote the part before the exponent to keep all the digits
      l <- gsub("^(.*)e", "'\\1'e", l)
      
      # turn the 'e+' into plotmath format
      l <- gsub("e", "%*%10^", l)
      # return this as an expression
      parse(text=l)
    }
    # Parameters table
    
    # Plotting the approach segment
    mytable <- data.frame(c(format(output$`Contact Point`,digits = 4),
                            format(output$`Zero force pt.`, digits = 4), 
                            format(output$`Young's Modulus`,digits = 4)),
                            row.names = c("C. Point:", "Zero F. Pt.:",
                                          "Young's Modulus:"))
    names(mytable) <- NULL
    # Theme for tableGrob
    mytheme <- ttheme_default(base_size = 10,core=list(fg_params=list(hjust=0, x=0.1)),
                              rowhead = list(fg_params = list(fontface = "bold")))
    # Base Plot
    applt <- plot(afmdata, segment = "approach") + 
      geom_vline(xintercept = output$`Contact Point`, lty = 2, lwd = .2)+
      geom_vline(xintercept = output$`Zero force pt.`, lty = 3, lwd = .3,
                 colour = "red")
    # Getting table position 
    xrange <- layer_scales(applt)$x$range$range[2] -
      layer_scales(applt)$x$range$range[1]
    xini <- layer_scales(applt)$x$range$range[1]
    xmin <- xini + 0.7*xrange
    xmax <- xini + 0.9*xrange
    
    yrange <- layer_scales(applt)$y$range$range[2] -
      layer_scales(applt)$y$range$range[1]
    yini <- layer_scales(applt)$y$range$range[1]
    ymin <- yini + 0.75*yrange
    ymax <- yini + 0.95*yrange
    # Add table to plot
    applt <- applt + annotation_custom(tableGrob(mytable, theme = mytheme), xmin = xmin , 
                              xmax = xmax, ymin = ymin, ymax = ymax) +
      ggtitle('Approach segment')
    
    # AÑADIR inset con la gráfica del ajuste de YM
    N <- nrow(afmdata$YoungModulus$fitdata)
    YMfitdata <- subset(afmdata$data, Segment == "approach" &
                          Z <= afmdata$Slopes$Z0Point, select = c(Z,ForceCorrected))
    YMfitdata <- YMfitdata[1:N,]
    applt <- applt + geom_point(data = YMfitdata, aes(x = Z, y = ForceCorrected),
                                colour = "red")
    inset <- ggplot(data = afmdata$YoungModulus$fitdata, 
                    aes(x = -Indentation, y = ForceCorrected)) + geom_point()
    inset <- inset + geom_line(aes(x = -Indentation, 
                                   y = afmdata$YoungModulus$fitYM$fitted.values), 
                               colour = "red", lwd = 1) + theme_bw()+
      scale_x_continuous(name = "Indentation",labels = fancy_scientific)+
      scale_y_continuous(name = "Force",labels = fancy_scientific)+
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 7))+
      theme(axis.text.y = element_text(size = 7))
    ins <- ggplotGrob(inset)
    applt + annotation_custom(grob = ins, 
                              xmin = xini + 0.4*xrange, 
                              xmax = xini + 1*xrange, 
                              ymax = yini + 0.75*yrange, 
                              ymin = yini + 0.1*yrange)
    if ("retract" %in% levels(afmdata$data$Segment)){
      
    }
  }
}