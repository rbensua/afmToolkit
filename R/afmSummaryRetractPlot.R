# Internal module of summary.afmdata(): plots the retract segment (when
# present) with the adhesion energy polygons and the energies table
# (when the adhesion analysis has been run).
afmSummaryRetractPlot <- function(afmdata) {
  Z <- ForceCorrected <- Segment <- NULL
  if (!("retract" %in% levels(afmdata$data$Segment))) {
    return(invisible(NULL))
  }
  retplt <- plot(afmdata, segment = "retract") + ggtitle("Retract segment")
  if ("AdhEner" %in% names(afmdata)) {
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
    xmin <- xini + 0.7 * xrange
    xmax <- xini + 0.9 * xrange
    yini <- layer_scales(retplt)$y$range$range[1]
    ymin <- yini + 0.3 * yrange
    ymax <- yini + 0.9 * yrange
    mytheme <- ttheme_default(base_size = 10,
                              core = list(fg_params = list(hjust = 0, x = 0.1)),
                              rowhead = list(fg_params = list(fontface = "bold")))
    
    afmdata$data %>% filter(Segment == "retract") %>%
      select(one_of("Z","ForceCorrected")) %>%
      filter(Z >= Z[afmdata$AdhEner$Points[1]] &
               Z <= Z[afmdata$AdhEner$Points[2]])  -> poly1
    afmdata$data %>% filter(Segment == "retract") %>%
      select(one_of("Z","ForceCorrected")) %>%
      filter(Z >= Z[afmdata$AdhEner$Points[2]] &
               Z <= Z[afmdata$AdhEner$Points[3]])  -> poly2
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
                 lty = 2, colour = "darkgrey") +
      annotation_custom(tableGrob(energies, theme = mytheme), xmin = xmin,
                        xmax = xmax, ymin = ymin, ymax = ymax)
  }
  print(retplt)
}
