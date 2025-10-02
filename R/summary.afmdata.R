#' @title Summary of an \code{afmdata} class object.
#'   
#'   
#' @description This function sumarizes the main features of an afmdata object and,
#' optionnaly plots all segments available with all parameters estimated.
#' 
#' @param object An object of \code{afmdata} class.
#' @param plt Logical variable. If TRUE plots all available segments with all available
#'   data.
#' @param ... Additional arguments (for compatibility with \code{summary})
#' @examples
#' \dontrun{path <- path.package("afmToolkit")
#' data <- afmReadJPK("force-save-JPK-3h.txt.gz", path = path)
#' data <- afmContactPoint(data, width = 20, mul1 = 1, mul2 = 10)
#' data <- afmDetachPoint(data, width = 20, mul1 = 2, mul2 = 30)
#' data <- afmBaselineCorrection(data)
#' data <- afmAdhesionEnergy(data, width = 20, mul = 10)
#' data <- afmZeroPointSlope(data, segment = "approach")
#' data <- afmIndentation(data)
#' data <- afmYoungModulus(data, thickness = 1e-7, params = list(alpha = 22),
#'                         silent = TRUE)
#' data <- afmExpDecay(data, nexp = 2, type = "CH")
#' summary(data)}                        
#' @import grid
#' @importFrom  gridExtra ttheme_default tableGrob
#' @import scales
#' @import dplyr
#' @importFrom graphics curve plot
#' @importFrom stats coefficients complete.cases
#' @importFrom utils read.table
#' @rdname summary
#' @export
#' @method summary afmdata

summary.afmdata <- function(object , plt = TRUE, ...){
  Indentation <- Time <- Force <- Segment <- Z <- ForceCorrected <- NULL
  afmdata <- object
  if (!is.afmdata(afmdata)){
    stop("Input should be an afmdata class object!")
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
    if(grepl("Z ~",afmdata$ExpFit$expdecayModel$call[2])){
      typeoffit <- "Constant Force"
    }else{
      typeoffit <- "Constant Height"
    }
    newnames <- c(names(output),"Type of Experiment")
    output <- cbind(output,typeoffit)
    names(output) <- newnames
    expfit <- summary(afmdata$ExpFit$expdecayModel)
    print(expfit$coefficients)
  }
  print(output)
  if (plt){
     
    fancy_scientific <- function(l) {
      # turn in to character string in scientific notation
      #print(l)
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
    applt <- plot(afmdata, segment = "approach") 
    if ("CP" %in% fields){
      applt<- applt+geom_vline(xintercept = output$`Contact Point`, lty = 2, lwd = .2)
    }
    if ("Slopes" %in% names(afmdata)){
      applt <- applt + geom_vline(xintercept = output$`Zero force pt.`, lty = 3, lwd = .3,
                 colour = "red")
    }
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
    
    if ("YoungModulus" %in% names(afmdata)){
    # add inset with plot of YM fit
    N <- nrow(afmdata$YoungModulus$fitdata)
    YMfitdata <- subset(afmdata$data, Segment == "approach" &
                          Z <= afmdata$Slopes$Z0Point, select = c(Z,ForceCorrected))
    YMfitdata <- YMfitdata[1:N,]
    applt <- applt + geom_point(data = YMfitdata, aes(x = Z, y = ForceCorrected),
                                colour = "red")
    inset <- ggplot(data = afmdata$YoungModulus$fitdata, 
                    aes(x = abs(Indentation), y = ForceCorrected)) + geom_point()
    inset <- inset + geom_line(aes(x = abs(Indentation), 
                                   y = afmdata$YoungModulus$fitYM$fitted.values), 
                               colour = "red", lwd = 1) + theme_bw()+
      scale_x_continuous(name = "Indentation") + #,labels = fancy_scientific)+
      scale_y_continuous(name = "Force") + #,labels = fancy_scientific)+
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 7))+
      theme(axis.text.y = element_text(size = 7))
    ins <- ggplotGrob(inset)
    applt <- applt + annotation_custom(grob = ins, 
                              xmin = xini + 0.5*xrange, 
                              xmax = xini + 1*xrange, 
                              ymax = yini + 0.75*yrange, 
                              ymin = yini + 0.2*yrange)
    }
    print(applt)
    if ("pause" %in% levels(afmdata$data$Segment)){
     ctc1 <- plot(afmdata, segment = "pause", vs = "Time") + xlab(NULL) + 
       scale_x_continuous(labels = NULL)+
       #scale_y_continuous(labels = fancy_scientific) +
       ggtitle("Pause segment")
     ctc2 <- ggplot(data = subset(afmdata$data, Segment == "pause"),
                    aes(x = Time, y = Z)) + geom_line() +
       #scale_y_continuous(labels = fancy_scientific)+
       theme_bw()
     if ("ExpFit" %in% names(afmdata)){
       times <- subset(afmdata$data, 
                       Segment == "pause")$Time[seq_along(afmdata$ExpFit$expdecayFit)]
       fitcoef <- summary(afmdata$ExpFit$expdecayModel)$coefficients[,1:2]
       fitcoef <- format(fitcoef, digits = 2)
       xrange <- layer_scales(ctc2)$x$range$range[2] -
         layer_scales(ctc2)$x$range$range[1]
       xini <- layer_scales(ctc2)$x$range$range[1]
       xmin <- xini + 0.7*xrange
       xmax <- xini + 0.9*xrange
       mytheme <- ttheme_default(base_size = 10,core=list(fg_params=list(hjust=0, x=0.1)),
                                 rowhead = list(fg_params = list(fontface = "bold")))
       if (typeoffit == "Constant Force"){
         yrange <- layer_scales(ctc2)$y$range$range[2] -
           layer_scales(ctc2)$y$range$range[1]
         yini <- layer_scales(ctc2)$y$range$range[1]
         ymin <- yini + 0.3*yrange
         ymax <- yini + 0.9*yrange
         ctc2 <- ctc2 + geom_line(data = data.frame(Time = times, 
                                                    Z = afmdata$ExpFit$expdecayFit),
                                  aes(x = Time, y = Z), col = "green", size = 1.5) +
           annotation_custom(tableGrob(fitcoef, theme = mytheme), xmin = xmin , 
                             xmax = xmax, ymin = ymin, ymax = ymax) 
       }else {
         yrange <- layer_scales(ctc1)$y$range$range[2] -
           layer_scales(ctc1)$y$range$range[1]
         yini <- layer_scales(ctc1)$y$range$range[1]
         ymin <- yini + 0.3*yrange
         ymax <- yini + 0.9*yrange
         ctc1 <- ctc1 + geom_line(data = data.frame(Time = times, 
                                                    Force = afmdata$ExpFit$expdecayFit),
                                  aes(x = Time, y = Force), col = "green", size = 1.5)+
           annotation_custom(tableGrob(fitcoef, theme = mytheme), xmin = xmin , 
                             xmax = xmax, ymin = ymin, ymax = ymax) 
       }
     }
       g1 <- ggplotGrob(ctc1)
       g2 <- ggplotGrob(ctc2)
       g <- rbind(g1,g2, size = "first")
       g$widths <- unit.pmax(g1$widths,g2$widths)
       g$layout[grepl("guide", g$layout$name),c("t","b")] <- c(1,nrow(g))
       grid.newpage()
       grid.draw(g)
    }
    if ("retract" %in% levels(afmdata$data$Segment)){
      retplt <- plot(afmdata, segment = "retract") + ggtitle("Retract segment")
      if ("AdhEner" %in% names(afmdata)){
        energies <- data.frame(abs(afmdata$AdhEner$Energies), row.names = NULL)
        colnames(energies) <- c("Adh. Energy", "Full detach Ener.", "Total Energy")
        energies <- format(energies, digits = 4)
        energies <- data.frame(t(energies))
        colnames(energies) <- NULL
        xrange <- layer_scales(retplt)$x$range$range[2] -
          layer_scales(retplt)$x$range$range[1]
        yrange <- layer_scales(retplt)$y$range$range[2] -
          layer_scales(retplt)$y$range$range[1]
        
        xini <- layer_scales(retplt)$x$range$range[1]
        xmin <- xini + 0.7*xrange
        xmax <- xini + 0.9*xrange
        yini <- layer_scales(retplt)$y$range$range[1]
        ymin <- yini + 0.3*yrange
        ymax <- yini + 0.9*yrange
        mytheme <- ttheme_default(base_size = 10,core=list(fg_params=list(hjust=0, x=0.1)),
                                  rowhead = list(fg_params = list(fontface = "bold")))
        
        
        
        afmdata$data %>% filter(Segment == "retract") %>% 
          select(one_of("Z","ForceCorrected")) %>% 
          filter(Z>=Z[afmdata$AdhEner$Points[1]] & 
                   Z<=Z[afmdata$AdhEner$Points[2]])  -> poly1
        afmdata$data %>% filter(Segment == "retract") %>% 
          select(one_of("Z","ForceCorrected")) %>% 
          filter(Z>=Z[afmdata$AdhEner$Points[2]] & 
                   Z<=Z[afmdata$AdhEner$Points[3]])  -> poly2
        poly1 <- bind_rows(data.frame(Z = poly1$Z[1], 
                                      ForceCorrected = 0), poly1)
        poly1 <- bind_rows(poly1, data.frame(Z = poly1$Z[nrow(poly1)], 
                                             ForceCorrected = 0))
        poly2 <- bind_rows(data.frame(Z = poly2$Z[1], 
                                      ForceCorrected = 0), poly2)
        poly2 <- bind_rows(poly2, data.frame(Z = poly2$Z[nrow(poly2)], 
                                             ForceCorrected = 0))
        retplt <- retplt + 
          geom_polygon(data = poly1, aes(x = Z, y = ForceCorrected), fill = "red") +
          geom_polygon(data = poly2, aes(x = Z, y = ForceCorrected), fill = "green") +
          geom_vline(xintercept = 
                       subset(afmdata$data, 
                              Segment == "retract")$Z[afmdata$AdhEner$Points],
                     lty = 2, colour = "darkgrey")+
          annotation_custom(tableGrob(energies, theme = mytheme), xmin = xmin , 
                            xmax = xmax, ymin = ymin, ymax = ymax) 
      }
      print(retplt)
    }
    
  }
  
     
}
