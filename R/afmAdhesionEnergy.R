afmAdhesionEnergy <-
  function(afmdata,
           width = 1,
           lagdiff = width,
           mul,
           mdj = NULL) {
    # Extract the retract segment from the afmdata
    if (is.afmexperiment(afmdata)) {
      AdhE <-
        lapply(afmdata, function(x)
          afmAdhesionEnergy(
            x,
            width = width,
            mul = mul,
            lagdiff = lagdiff,
            mdj = mdj
          ))
      afmexperiment <-
        mapply(
          afmdata,
          AdhE,
          FUN = function(x, y)
            append.afmdata(x, y, name = "AdhEner"),
          SIMPLIFY = FALSE
        )
      return(afmexperiment(afmexperiment))
    } else if (is.afmdata(afmdata)) {
      data <- subset(afmdata$data, Segment == "retract")
      Z <- data$Z
      Force <- data$ForceCorrected
      
      
      imin <- min(which(Force < 0))
      idx <- seq_along(Force)
      imax <- min(which(Force > 0 & idx > which.min(Force)))
      
      
      
      
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
      
      #Noise at the begining of the curve
      noiseEndDelta <-
        sd(delta[width:(width + as.integer(0.1 * n))], na.rm =
             TRUE)
      # Noise at the end of the curve
      noiseBeginDelta <-
        sd(delta[(n - (width + as.integer(0.1 * n))):(n - width)],
           na.rm = TRUE)
      # Tolerances in Delta
      
      tol2Begin <- mul * noiseBeginDelta
      tol2End <- mul * noiseEndDelta
      
      # Jumps in Delta
      jumpsDelta <- list()
      
      # Indices above the threshold (in the slope)
      idxGrTol2Begin <- which(abs(delta) > tol2Begin)
      idxGrTol2Begin <- idxGrTol2Begin[idxGrTol2Begin <= imax + width &
                                         idxGrTol2Begin >= which.min(Force)]
      # Compute the number of jumps at the beginning
      
      Njumps <- sum(diff(idxGrTol2Begin) > mdj) + 1
      idxgr <- which(diff(idxGrTol2Begin) > mdj)
      
      gr.begin <- c(idxGrTol2Begin[1], idxGrTol2Begin[idxgr + 1])
      gr.end <- c(idxGrTol2Begin[idxgr], tail(idxGrTol2Begin, n = 1))
      
      
      for (k in 1:Njumps) {
        Jump <- gr.begin[k]:gr.end[k]
        jumpsDelta <- append(jumpsDelta, list(Jump))
      }
      
      
      
      # Magnitudes of the jumps (in 2 ways)
      
      magnitude1 <- numeric() # Computing the range
      magnitude1Corrected <- numeric()
      magnitude2 <- numeric() # Integrating bRoll
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
        windowJump <- round(pmedio - radius2):round(pmedio + radius2)
        magnitude1[i] <-
          diff(range(Force[(pmedio - radius1):(pmedio + radius1)])) -
          2 * (sqrt(2 * log(width)) - 0.588) * noiseForce
        
        magnitude1Corrected[i] <-
          diff(range(Force[(pmedio - radius1):(pmedio + radius1)] -
                       (bRoll[round(pmedio + radius2)] +
                          bRoll[round(pmedio - radius2)]) / 2 *
                       Z[(pmedio - radius1):(pmedio + radius1)])) -
          2 * (sqrt(2 * log(width)) - 0.588) * noiseForce
        
        zWJ <- Z[windowJump]
        bWJCorrected <-
          bRoll[windowJump] - (bRoll[round(pmedio + radius2)] +
                                 bRoll[round(pmedio - radius2)]) / 2
        magnitude2[i] <- trapz(zWJ, bRoll[windowJump])
        magnitude2Corrected[i] <- trapz(zWJ, bWJCorrected)
      }
      magnitudeDF <- data.frame(
        pmid = pmid,
        magnitude1 = magnitude1,
        magnitude1Corrected = magnitude1Corrected,
        magnitude2 = magnitude2,
        magnitude2Corrected = magnitude2Corrected
      )
      
      bigJump <- which.max(magnitudeDF$magnitude2Corrected)
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
      return(
        list(
          Magnitudes = magnitudeDF,
          Points = c(imin, leftPoint, imax),
          delta = delta,
          noiseB = noiseBeginDelta,
          noiseE = noiseEndDelta,
          EnerFun = data.frame(Zrel = Zadh - max(Zadh), EnerFun = EnerFun),
          Energies = data.frame(
            E1adh = E1adh,
            E2adh = E2adh,
            Weight = Weight
          ),
          jumpsDelta = jumpsDelta
        )
      )
    } else{
      stop("No afmdata or afmexperiment class input provided.")
    }
  }
