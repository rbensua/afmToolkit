<<<<<<< HEAD
isPeak2 <- function (f, SoN = 2, span = 81, sm.span = 11, plot = FALSE,
          add = FALSE, zerothrsh = 2, area.w = 0.003, ratio = 0.2,
          ...)
{
  parea <- function(f, pt, eps = area.w) {
    x <- f[, 1]
    sum(f[x <= pt * (1 + eps) & x >= pt * (1 - eps), 2])
  }
  n <- dim(f)[1]
  mz <- f[, 1]
  lo <- lnn(f[, 2], span = span, sm.span = sm.span)
  sm <- lo$fitted
  ispeak <- peaks(sm, span)
  sigmas <- lo$s
  peak <- ispeak & sm > zerothrsh & sm > SoN * sigmas
  area <- sapply(mz[peak], parea, f = cbind(mz, sm)[!peak,
                                                    ])
  As <- rep(NA, length(peak))
  As[peak] <- area
  peak[peak] <- peak[peak] & area/max(area) > ratio
  pks <- data.frame(peak = peak, smooth = sm, mz = mz, sigmas = sigmas,
                    area = As)
  if (plot) {
    if (add) {
      lines(mz, sm, col = "cyan")
      points(mz[peak], sm[peak], col = "orange")
    }
    else specZoom(pks, ...)
  }
  return(pks)
}
=======
isPeak2 <- function (f, SoN = 2, span = 81, sm.span = 11, plot = FALSE,
          add = FALSE, zerothrsh = 2, area.w = 0.003, ratio = 0.2,
          ...)
{
  parea <- function(f, pt, eps = area.w) {
    x <- f[, 1]
    sum(f[x <= pt * (1 + eps) & x >= pt * (1 - eps), 2])
  }
  n <- dim(f)[1]
  mz <- f[, 1]
  lo <- lnn(f[, 2], span = span, sm.span = sm.span)
  sm <- lo$fitted
  ispeak <- peaks(sm, span)
  sigmas <- lo$s
  peak <- ispeak & sm > zerothrsh & sm > SoN * sigmas
  area <- sapply(mz[peak], parea, f = cbind(mz, sm)[!peak,
                                                    ])
  As <- rep(NA, length(peak))
  As[peak] <- area
  peak[peak] <- peak[peak] & area/max(area) > ratio
  pks <- data.frame(peak = peak, smooth = sm, mz = mz, sigmas = sigmas,
                    area = As)
  if (plot) {
    if (add) {
      lines(mz, sm, col = "cyan")
      points(mz[peak], sm[peak], col = "orange")
    }
    else specZoom(pks, ...)
  }
  return(pks)
}
>>>>>>> 74325c08fcfc012ad5d5e5ed8c461a93115269f1
