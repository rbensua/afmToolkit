# Internal module of summary.afmdata(): plots the approach segment with
# the parameter table, the contact point / zero force point vertical
# lines and (when present) the Young's modulus fit inset.
afmSummaryApproachPlot <- function(afmdata, output) {
  Indentation <- Z <- ForceCorrected <- Segment <-  NULL
  # Parameters table: only the parameters actually present are shown
  # (a missing field would leave empty entries, which used to break the
  # construction of the table)
  table_entries <- c(
    if (!is.null(output$`Contact Point`))
      c("C. Point:" = format(output$`Contact Point`, digits = 4)),
    if (!is.null(output$`Zero force pt.`))
      c("Zero F. Pt.:" = format(output$`Zero force pt.`, digits = 4)),
    if (!is.null(output$`Young's Modulus`))
      c("Young's Modulus:" = format(output$`Young's Modulus`, digits = 4))
  )
  
  # Theme for tableGrob
  mytheme <- ttheme_default(base_size = 10,
                            core = list(fg_params = list(hjust = 0, x = 0.1)),
                            rowhead = list(fg_params = list(fontface = "bold")))
  # Base Plot
  applt <- plot(afmdata, segment = "approach")
  if ("CP" %in% names(afmdata)) {
    applt <- applt + geom_vline(xintercept = output$`Contact Point`, lty = 2, lwd = .2)
  }
  if ("Slopes" %in% names(afmdata)) {
    applt <- applt + geom_vline(xintercept = output$`Zero force pt.`, lty = 3, lwd = .3,
                                colour = "red")
  }
  # Getting table position
  xrange <- layer_scales(applt)$x$range$range[2] -
    layer_scales(applt)$x$range$range[1]
  xini <- layer_scales(applt)$x$range$range[1]
  xmin <- xini + 0.7 * xrange
  xmax <- xini + 0.9 * xrange
  
  yrange <- layer_scales(applt)$y$range$range[2] -
    layer_scales(applt)$y$range$range[1]
  yini <- layer_scales(applt)$y$range$range[1]
  ymin <- yini + 0.75 * yrange
  ymax <- yini + 0.95 * yrange
  # Add table to plot (only if there is something to show in it)
  if (length(table_entries) > 0) {
    mytable <- data.frame(table_entries, row.names = names(table_entries))
    names(mytable) <- NULL
    applt <- applt +
      annotation_custom(tableGrob(mytable, theme = mytheme), xmin = xmin,
                        xmax = xmax, ymin = ymin, ymax = ymax)
  }
  applt <- applt + ggtitle('Approach segment')
  
  if ("YoungModulus" %in% names(afmdata)) {
    # add inset with plot of YM fit
    N <- nrow(afmdata$YoungModulus$fitdata)
    YMfitdata <- subset(afmdata$data, Segment == "approach" &
                          Z <= afmdata$Slopes$Z0Point,
                        select = c(Z, ForceCorrected))
    YMfitdata <- YMfitdata[1:N,]
    applt <- applt + geom_point(data = YMfitdata,
                                aes(x = Z, y = ForceCorrected),
                                colour = "red")
    inset <- ggplot(data = afmdata$YoungModulus$fitdata,
                    aes(x = abs(Indentation), y = ForceCorrected)) + geom_point()
    inset <- inset +
      geom_line(aes(x = abs(Indentation),
                    y = afmdata$YoungModulus$fitYM$fitted.values),
                colour = "red", lwd = 1) + theme_bw() +
      scale_x_continuous(name = "Indentation") +
      scale_y_continuous(name = "Force") +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 7)) +
      theme(axis.text.y = element_text(size = 7))
    ins <- ggplotGrob(inset)
    applt <- applt + annotation_custom(grob = ins,
                                       xmin = xini + 0.5 * xrange,
                                       xmax = xini + 1 * xrange,
                                       ymax = yini + 0.75 * yrange,
                                       ymin = yini + 0.2 * yrange)
  }
  print(applt)
}
