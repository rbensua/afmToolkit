windowed_scalogram_int <-
  function(signal, scales, windowrad, wname = c("MORLET","HAAR")) {
    wname <- toupper(wname)
    
    
    ns <- length(scales)
    nt <- length(signal)
    t_ini <- 1 + windowrad
    t_fin <- nt - windowrad
    ntw <- t_fin - t_ini + 1
    
    
    coefs <- cwt_optimized(scales, signal, wname)
    abscoefs2 <- abs(coefs) ^ 2
    
    if (wname == "MORLET") {
      waverad <- 3
      
      wr <- waverad * scales
      wrmax <- ceiling(wr[length(wr)])
      
      i_ini <- max(1, 1 + wrmax - 2 * windowrad)
      i_fin <- min(ntw, nt - wrmax)
      
      wsc <- matrix(0, nrow = i_fin - i_ini + 1, ncol = ns)
      
      for (i in 1:(i_fin - i_ini + 1)) {
        for (j in 1:ns) {
          wt_ini_j <- max(i + i_ini - 1, 1 + ceiling(wr[j]))
          
          wt_fin_j <-
            min(i + i_ini - 1 + 2 * windowrad, nt - ceiling(wr[j]))
          wsc[i, j] <- sum(abscoefs2[j, wt_ini_j:wt_fin_j])
          # Normalization
          wsc[i, j] <-
            wsc[i, j] * (2 * windowrad + 1) / (wt_fin_j - wt_ini_j +
                                                 1)
        }
      }
    } else if (wname == "HAAR") {
      waverad <- 3 / 4
      
      wr <- waverad * scales
      wrmax <- ceiling(wr[length(wr)])
      
      i_ini <- max(1, 1 + wrmax - 2 * windowrad)
      i_fin <- min(ntw, nt - wrmax)
      wsc <- matrix(0, nrow = i_fin - i_ini + 1, ncol = ns)
      for (i in 1:(i_fin - i_ini + 1)) {
        for (j in 1:ns) {
          wt_ini_j <- max(i + i_ini - 1, 1 + ceiling(wr[j]))
          
          wt_fin_j <-
            min(i + i_ini - 1 + 2 * windowrad, nt - ceiling(wr[j]))
          wsc[i, j] <- sum(abscoefs2[j, wt_ini_j:wt_fin_j])
          # Normalization
          wsc[i, j] <-
            wsc[i, j] * (2 * windowrad + 1) / (wt_fin_j - wt_ini_j +
                                                 1)
        }
      }
    } else{
      stop('At this time only Haar and Morlet are available...')
    }
    wsc <- sqrt(wsc)
    tt_ini <-  i_ini + windowrad
    tt_fin <- i_fin + windowrad
    return(list(
      wsc = wsc,
      tt_ini = tt_ini,
      tt_fin = tt_fin
    ))
  }