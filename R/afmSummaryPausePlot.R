# Internal module of summary.afmdata(): plots the pause segment (when
# present) as a two-panel figure (Force vs Time on top, Z vs Time below)
# with the relaxation fit overlay and its coefficients table.
# typeoffit ("Constant Force" / "Constant Height") comes from
# afmSummaryTable().
afmSummaryPausePlot <- function(afmdata, typeoffit = NULL) {
  Time <- Z <- Segment <- NULL
  if (!("pause" %in% levels(afmdata$data$Segment))) {
    return(invisible(NULL))
  }
  ctc1 <- plot(afmdata, segment = "pause", vs = "Time") + xlab(NULL) +
    scale_x_continuous(labels = NULL) +
    ggtitle("Pause segment")
  ctc2 <- ggplot(data = subset(afmdata$data, Segment == "pause"),
                 aes(x = Time, y = Z)) + geom_line() +
    theme_bw()
  if ("Relax" %in% names(afmdata)) {
    times <- subset(afmdata$data,
                    Segment == "pause")$Time[seq_along(afmdata$Relax$relaxFit)]
    fitcoef <- summary(afmdata$Relax$relaxModel)$coefficients[,1:2]
    fitcoef <- format(fitcoef, digits = 2)
    xrange <- layer_scales(ctc2)$x$range$range[2] -
      layer_scales(ctc2)$x$range$range[1]
    xini <- layer_scales(ctc2)$x$range$range[1]
    xmin <- xini + 0.7 * xrange
    xmax <- xini + 0.9 * xrange
    mytheme <- ttheme_default(base_size = 10,
                              core = list(fg_params = list(hjust = 0, x = 0.1)),
                              rowhead = list(fg_params = list(fontface = "bold")))
    if (typeoffit == "Constant Force") {
      yrange <- layer_scales(ctc2)$y$range$range[2] -
        layer_scales(ctc2)$y$range$range[1]
      yini <- layer_scales(ctc2)$y$range$range[1]
      ymin <- yini + 0.3 * yrange
      ymax <- yini + 0.9 * yrange
      ctc2 <- ctc2 +
        geom_line(data = data.frame(Time = times,
                                    Z = afmdata$Relax$relaxFit),
                  aes(x = Time, y = Z), col = "green", linewidth = 1.1) +
        annotation_custom(tableGrob(fitcoef, theme = mytheme), xmin = xmin,
                          xmax = xmax, ymin = ymin, ymax = ymax)
    } else {
      yrange <- layer_scales(ctc1)$y$range$range[2] -
        layer_scales(ctc1)$y$range$range[1]
      yini <- layer_scales(ctc1)$y$range$range[1]
      ymin <- yini + 0.3 * yrange
      ymax <- yini + 0.9 * yrange
      ctc1 <- ctc1 +
        geom_line(data = data.frame(Time = times,
                                    Force = afmdata$Relax$relaxFit),
                  aes(x = Time, y = Force), col = "green", linewidth = 1.1) +
        annotation_custom(tableGrob(fitcoef, theme = mytheme), xmin = xmin,
                          xmax = xmax, ymin = ymin, ymax = ymax)
    }
  }
  g1 <- ggplotGrob(ctc1)
  g2 <- ggplotGrob(ctc2)
  g <- rbind(g1, g2, size = "first")
  g$widths <- unit.pmax(g1$widths, g2$widths)
  g$layout[grepl("guide", g$layout$name), c("t","b")] <- c(1, nrow(g))
  grid.newpage()
  grid.draw(g)
}
