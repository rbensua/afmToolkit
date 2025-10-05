#' @title Adhesion Energy
#'
#' @description
#' Finds the adhesion and the full detach energies from the retract segment of the AFM 
#' F-d curve. 
#' @usage afmAdhesionEnergy(afmdata, width = 1, lagdiff = width, mul, mdj = NULL)
#' @param afmdata An afmdata or afmexperiment class variables. Baseline correction should
#'  have been done already.
#' @param width Width of the window for the local regression (in vector position units)
#' @param lagdiff Lag for estimating the differences in Delta (or slopes) signal. 
#' By default it takes the same value as the window with.
#' @param mul Multiplier for the calculating the threshold inthe estimation of jumps 
#' and peaks in the Delta signal
#' @param mdj Minimum distance between jumps. If none is given then it will be set equal 
#' to \code{width}
#' @return An \code{afmdata} class variable which will consist on the original 
#' input \code{afmdata} variable plus a new list named \code{AdhEner} with the 
#' following fields:
#'
#' \code{Points} Array containing the indices of the retract segment where the adhesion 
#' begins, the unbinding event takes place and the adhesion ends.
#'
#' \code{Energies} Data frame with three columns: \code{E1adh}, \code{E2adh} and 
#' \code{Etotal}, being the first one the energy from the begining of the adhesion until
#' the unbinding event, then second one the energy from the unbinding event until the full
#' detachment of the tip, and the third one, the sum of them.
#' @examples
#' path <- path.package("afmToolkit")
#' data <- afmReadJPK("force-save-JPK-3h.txt.gz", path = path)
#' data <- afmContactPoint(data, width = 20, mul1 = 1, mul2 = 10)
#' data <- afmDetachPoint(data, width = 20, mul1 = 2, mul2 = 30)
#' data <- afmBaselineCorrection(data)
#' data <- afmAdhesionEnergy(data, width = 20, mul = 10)
#' str(data$AdhEner)
#' @importFrom stats sd
#' @importFrom utils tail
#' @export
afmAdhesionEnergy <-
  function(afmdata,
           width = 1,
           lagdiff = width,
           mul,
           mdj = NULL) {
    # Extract the retract segment from the afmdata
    Segment <- NULL
    if (is.afmexperiment(afmdata)) {
      afmexperiment <-
        lapply(afmdata, function(x) {
          if (!is.null(x$params$curvename)) {
            print(paste("Processing curve: ", x$params$curvename), sep = " ")
          }
          afmAdhesionEnergy(
            x,
            width = width,
            mul = mul,
            lagdiff = lagdiff,
            mdj = mdj
          )
        })
      return(afmexperiment(afmexperiment))
    } else if (is.afmdata(afmdata)) {
      data <- subset(afmdata$data, Segment == "retract")
      Z <- data$Z
      Force <- data$ForceCorrected
      
      
      imin <- min(which(Force < 0))
      idx <- seq_along(Force)
      iFmin <- which.min(Force)
      Fmin <- min(Force)
      ZFmin <- Z[iFmin]
      imax <- min(which(Force > 0 & idx > iFmin))
      
      
      
      
      n = length(Z)
      if (is.null(mdj)) {
        mdj = width
      }
      
      delta <- array(0, dim = c(n, 1))
      app <- matrix(c(rep(1, n), Z, Force), nrow = n, ncol = 3)
      bRoll <- windowedFit(app, width)
      dif1 <- diff(bRoll, lag = lagdiff)
      delta <- dif1
      delta <- c(rep(0, width + lagdiff), delta, rep(0, width))
      bRoll <- c(rep(0, width), bRoll, rep(0, width))
      
      #Noise at the contact part of the curve
#      noiseContactDelta <-
#        sd(delta[width:(width + as.integer(0.1 * n))], na.rm =
#             TRUE)

            # Noise at the non-contact part of the curve
      noiseFreeDelta <-
        sd(delta[(n - (width + as.integer(0.1 * n))):(n - width)],
           na.rm = TRUE)
      # Tolerances in Delta
      
      tolFree <- mul * noiseFreeDelta
 #     tolContact <- mul * noiseContactDelta
      
      # Jumps in Delta
      jumpsDelta <- list()
      
      # Indices above the threshold (in the slope)
      idxGrTolFree <- which(abs(delta) > tolFree)
      idxGrTolFree <- idxGrTolFree[idxGrTolFree <= imax + width &
                                       idxGrTolFree >= which.min(Force)]
      # Compute the number of jumps at the beginning
      
      Njumps <- sum(diff(idxGrTolFree) > mdj) + 1
      idxgr <- which(diff(idxGrTolFree) > mdj)
      
      gr.begin <- c(idxGrTolFree[1], idxGrTolFree[idxgr + 1])
      gr.end <-
        c(idxGrTolFree[idxgr], tail(idxGrTolFree, n = 1))
      
      
      for (k in 1:Njumps) {
        Jump <- gr.begin[k]:gr.end[k]
        jumpsDelta <- append(jumpsDelta, list(Jump))
      }
      
      
      
      # Magnitudes of the jumps
      
      
      magnitude2Corrected <- numeric()
      noiseForce <-
        sd(Force[(n - (width + as.integer(0.1 * n))):(n - width)],
           na.rm = TRUE)
      pmid <- numeric()
      
      for (i in seq_along(jumpsDelta)) {
        r1 <- 1.25 * width
        r2 <- 0.5 * (length(jumpsDelta[[i]]) - lagdiff)
        radius1 <- max(r2, width)
        radius2 <- max(r1, r2)
        pmedio <-
          as.integer(jumpsDelta[[i]][1] + tail(jumpsDelta[[i]], n = 1)) /
          2 - lagdiff / 2
        pmid[i] <- pmedio
        windowJump <-
          round(pmedio - radius2):round(pmedio + radius2)
        zWJ <- Z[windowJump]
        bWJCorrected <-
          bRoll[windowJump] - (bRoll[round(pmedio + radius2)] +
                                 bRoll[round(pmedio - radius2)]) / 2
        magnitude2Corrected[i] <- trapz(zWJ, bWJCorrected)
      }
      magnitudeDF <- data.frame(pmid = pmid,
                                magnitude = magnitude2Corrected)
      
      bigJump <- which.max(magnitudeDF$magnitude)
      pmidbigJump <- pmid[bigJump]
      leftPoint <-
        pmidbigJump + max(0.5 * (length(jumpsDelta[[bigJump]]) - lagdiff), width)
      leftPoint <- min(leftPoint, imax)
      
      
      #Adhesion Energy
      Z1adh <- Z[imin:leftPoint]
      F1adh <- Force[imin:leftPoint]
      Z2adh <- Z[leftPoint:imax]
      F2adh <- Force[leftPoint:imax]
      E1adh <- trapz(Z1adh, F1adh)
      E2adh <- trapz(Z2adh, F2adh)
      if (is.na((E2adh))) {
        E2adh <- 0
      }
      Zadh <- Z[seq(from = imax, to = imin, by = -1)]
      Fadh <- Force[seq(from = imax, to = imin, by = -1)]
      Eadh <- numeric()
      for (i in seq_along(Zadh)) {
        Eadh[i] <- trapz(Zadh[1:i], Fadh[1:i])
      }
      Weight <- Eadh[length(Eadh)]
      EnerFun <- Eadh / Eadh[length(Eadh)]
      EnerFun[1] <- 0
      AdhEner <-  list(
        Points = c(imin, leftPoint, imax),
        Energies = data.frame(
          E1adh = abs(E1adh),
          E2adh = abs(E2adh),
          Etotal = Weight
        ),
        minForce = data.frame(
          Fmin = Fmin, 
          ZFmin = ZFmin
        )
      )
      return(append.afmdata(afmdata, AdhEner))
    } else{
      stop("No afmdata or afmexperiment class input provided.")
    }
  }
