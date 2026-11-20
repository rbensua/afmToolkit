#' @title Jump detector
#'
#' @description Jump detector for AFM curves. It locates abrupt jumps (tethers,
#' protein unfolding events, etc.) in an equally spaced signal, typically the
#' force of the approach or retract segment of a force-distance curve. Jumps
#' are found through the peaks of the Haar continuous wavelet transform of
#' the signal: positions with a coefficient energy above \code{mean + tol * sd}
#' are marked as significant, and consecutive significant positions closer
#' than \code{jumpsep} are grouped into a single jump, located at the local
#' maximum of the coefficient energy.
#'
#' @param signal An array with the values of the forces along time
#'   (equispaced).
#' @param jumpsep Minimum separation (in signal indices) between jumps.
#' @param jumpdur Maximum duration (in signal indices) of jumps. It is used to
#'   compute the maximum scale in the coefficients of the wavelet transform
#'   and as the half width of the window where the amplitude of each jump
#'   is measured.
#' @param tol Tolerance for the wavelet coefficients series to detect jumps,
#'   multiplied by the standard deviation of the wavelet coefficients series.
#' @param makefigure Logical. If TRUE, a two-panel figure is drawn: the
#'   wavelet coefficients with the detection threshold, and the signal with
#'   the amplitude of the detected jumps.
#'
#' @return A list with the following components:
#'   \itemize{
#'     \item \code{jump_idx}: indices of the signal where the jumps were
#'       detected (\code{NULL} if no jump was found).
#'     \item \code{Fmin}: minimum of the signal around each jump.
#'     \item \code{Fmax}: maximum of the signal around each jump.
#'   }
#'
#' @examples
#' # Synthetic signal with three jumps of similar magnitude
#' signal <- c(rep(0, 200), rep(1, 150), rep(0, 150), rep(1, 100))
#' jd <- jumpdetector(signal)
#' jd$jump_idx
#' jd$Fmax - jd$Fmin # amplitude of each jump
#'
#' @importFrom wavScalogram cwt_wst
#' @importFrom graphics abline lines par plot
#' @importFrom stats sd
#' @export
jumpdetector <- function(signal,
                         jumpsep = 40,
                         jumpdur = 4,
                         tol = 2,
                         makefigure = FALSE) {
  if (!is.numeric(signal)) {
    stop("'signal' must be a numeric vector.")
  }
  if (length(signal) == 0L || anyNA(signal)) {
    stop("'signal' must not be empty nor contain missing values.")
  }
  if (jumpsep < 1) {
    stop("'jumpsep' must be >= 1 (in signal indices).")
  }
  jumpdur <- floor(jumpdur)
  if (jumpdur < 1) {
    stop("'jumpdur' must be >= 1 (in signal indices).")
  }

  nt <- length(signal)
  minscale <- 2
  maxscale <- max(minscale, floor(2 * jumpdur))
  # The wavelet transform is computed with the power-of-two scales built
  # by cwt_wst() up to maxscale, so its border effects extend about
  # maxscale indices at each end of the signal. Those regions are
  # excluded, and the interior must remain non-empty for the detection
  # threshold to be meaningful.
  if (nt < 2 * maxscale + 2) {
    stop("signal is too short: it needs at least ", 2 * maxscale + 2,
         " points for jumpdur = ", jumpdur, ".")
  }

  kk <- cwt_wst(signal,
                scales = c(minscale, maxscale, 1),
                wname = "HAAR",
                makefigure = FALSE)
  wy <- rowSums(abs(kk$coefs))
  wy[1:maxscale] <- 0 # Border effects
  wy[(nt - maxscale):nt] <- 0
  nzwy <- wy[(maxscale + 1):(nt - maxscale - 1)]
  abstol <- mean(nzwy) + sd(nzwy) * tol

  jump_idx <- NULL
  Fmin <- NULL
  Fmax <- NULL

  signif_idx <- which(wy > abstol)
  nsi <- length(signif_idx)
  if (nsi > 0) {
    # Group the significant indices: a gap larger than jumpsep ends the
    # group. Each group is one jump, located at the local maximum of
    # the coefficient energy within the group.
    kkjump <- c(which(diff(signif_idx) > jumpsep), nsi)
    njump <- length(kkjump)
    jump_idx <- integer(njump)
    jend <- 0L
    for (i in 1:njump) {
      jini <- jend + 1L
      jend <- kkjump[i]
      group <- signif_idx[jini:jend]
      jump_idx[i] <- group[which.max(wy[group])]
    }

    # Amplitude of each jump: range of the signal in a window of
    # +/- jumpdur indices around the jump.
    Fmin <- numeric(njump)
    Fmax <- numeric(njump)
    for (i in 1:njump) {
      jini <- max(1, jump_idx[i] - jumpdur)
      jfin <- min(nt, jump_idx[i] + jumpdur)
      Fmin[i] <- min(signal[jini:jfin])
      Fmax[i] <- max(signal[jini:jfin])
    }
  }

  if (makefigure) {
    # A single two-panel figure (the standalone version drew two
    # separate devices) with the user's par() restored on exit.
    op <- par(mfrow = c(1, 2), no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    mnzwy <- mean(nzwy)
    plot(wy, type = "l",
         xlab = "Index", ylab = "Wavelet coefficients",
         main = "Wavelet coefficients and detection threshold")
    abline(h = mnzwy, lty = "dotted")
    abline(h = mnzwy + sd(nzwy), lty = "dotdash")
    abline(h = abstol)
    plot(signal, type = "l",
         xlab = "Index", ylab = "Signal",
         main = "Signal and detected jumps")
    if (!is.null(jump_idx)) {
      for (i in 1:njump) {
        lines(c(jump_idx[i], jump_idx[i]), c(Fmin[i], Fmax[i]),
              col = "red", lwd = 2)
      }
    }
  }

  list(jump_idx = jump_idx, Fmin = Fmin, Fmax = Fmax)
}
