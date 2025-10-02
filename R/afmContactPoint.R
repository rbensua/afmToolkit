#' @title Contact point
#'
#' @description
#' Find the contact point in for the Force-Distance curve
#' following the local regression and two thresholds methods described 
#' in Microscopy Research and Technique 2013 (see reference).
#' 
#' @param afmdata A Force-Distance curve with the afmdata structure. It should be a list
#'   with at least the 'data' field with a data frame of at least 4 columns.
#' @param width Width of the window for the local regression (in vector position units)
#' @param mul1 First multiplier for the first alarm threshold
#' @param mul2 Second multiplier for the second alarm threshold
#' @param lagdiff Lag for estimating the differences in Delta (or slopes) signal. By
#'   default it takes the same value as the window with.
#' @param Delta Logical. If TRUE, then the statistic for determining the contact point is
#'   the differences between two consecutive values of the slope of the local regression
#'   line. If FALSE then the slope itself is used.
#' @param loessSmooth Logical If TRUE, a loess smoothing (via loess.smooth()) is done
#'   prior to the determination of the contact point. The span of the  smoothing is 0.05
#'   (5%), the degree is 2 and the number of points equals the number of points in the
#'   approach segment.
#' @param silent Logical. If TRUE supresses a message showing which curve is being processing
#'    (use it when processing a high number of curves). Defaults to FALSE.
#' @return An \code{afmdata} class variable which will consist on the original input
#'   \code{afmdata} variable plus a new list named \code{CP} with the following fields:
#'
#' \code{CP} The contact point value.
#'
#' \code{iCP} The position in the array for the contact point value.
#'
#' \code{delta} The delta signal.
#'
#' \code{noise} The noise of the delta signal
#' @examples
#' path <- path.package("afmToolkit")
#' data <- afmReadJPK("force-save-JPK-3h.txt.gz", path = path)
#' width <- 20
#' mul1 <- 1
#' mul2 <- 10
#' data <- afmContactPoint(data, width = width, mul1 = mul1, mul2 = mul2)
#' \dontrun{
#' plot(data, segment = "approach") + geom_vline(xintercept = data$CP$CP, lty = 2)
#' }
#' @references 
#'  Benitez R., Moreno-Flores S., Bolos V. J. and Toca-Herrera J.L. (2013). "A 
#'  new automatic contact point detection algorithm for AFM force curves". 
#'  Microscopy research and technique, \strong{76} (8), pp. 870-876.
#' @seealso \code{\link{afmDetachPoint}}
#' @importFrom stats loess.smooth
#' @export



afmContactPoint <-
  function(afmdata,
           width = 1,
           mul1,
           mul2,
           lagdiff = width,
           Delta = TRUE,
           loessSmooth = FALSE, 
           silent = FALSE) {
    
    Segment <- NULL
# Checking data inputs ####
    if(is.afmexperiment(afmdata)) {
      #* Case afmexperiment ------
      afmexperiment <- lapply(afmdata, function(x) {
        if (!is.null(x$params$curvename) & !silent) {
          print(paste("Processing curve: ", x$params$curvename), sep = " ")
        }
        afmContactPoint(
          x,
          width = width,
          mul1 = mul1,
          mul2 = mul2,
          lagdiff = lagdiff,
          Delta = Delta,
          loessSmooth = loessSmooth
        )
      })
      return(afmexperiment(afmexperiment))
    }else if(is.afmdata(afmdata)){
      #* Case afmdata 1 curve only ----
      {
      # If afmdata is a multi-indentation experiment, we obtain the contact 
      # point with the first approach segment.
      approach_string <- ifelse(is.afmmulti(afmdata),"approach1","approach")
      data.approach <- subset(afmdata$data, Segment == approach_string)
      n <- nrow(data.approach)
      direction <- data.approach$Z[n] - data.approach$Z[1]
      if (loessSmooth) {
        #* Smoothing the data ----
        data.approach.smoothed <-
          loess.smooth(
            data.approach$Z,
            data.approach$Force,
            span = 0.05,
            degree = 2,
            evaluation = n
          )
        Z <- data.approach.smoothed$x
        Force <- data.approach.smoothed$y
        if (direction < 0) {
          Z <- rev(Z)
          Force <- rev(Force)
        }
      } else{
        #* No smoothing ----
        Z <- data.approach$Z
        Force <- data.approach$Force
      }
      
      #* Creating the arrays: bRoll -> slopes in a rolling windows and delta -> change in slopes. -----
      delta <- array(0, dim = c(n, 1))
      imax <- n - width
      app <- matrix(c(rep(1, n), Z, Force), nrow = n, ncol = 3)
      bRoll <- windowedFit(app, width)
      delta <- diff(bRoll, lag = lagdiff)
      delta <- c(rep(0, width + lagdiff), delta, rep(0, width))
      if (!Delta) {
        delta <- c(rep(0, width), bRoll, rep(0, width))
      }
      
      noise <-
        sd(delta[(width + as.integer(n / 3)):(width + as.integer(n / 3) + 
                                                as.integer(0.1 * n))],
           na.rm = TRUE)
      tol1 <- mul1 * noise
      tol2 <- mul2 * noise
      
      if (tol2 > max(abs(delta))) {
        tol2 <- max(abs(delta)) - 0.05 * diff(range(abs(delta)))
      }
      
      idxGrTol2 <- which(abs(delta) > tol2)
      idxSmTol1 <- which(abs(delta) < tol1)

      j <- max(idxSmTol1[idxSmTol1 < min(idxGrTol2)])#+1
      
      if ((j > 1) & (delta[j] != 0)) {
        eps <-
          (tol1 - abs(delta[j])) / abs(delta[j + 1] - delta[j]) 
      }
      else{
        eps <- 0
      }
      i_contact = min(j + width, n - 1)
      z_contact = Z[i_contact]
      z_contact = z_contact + eps * (Z[i_contact + 1] - z_contact)
      
      CP <- list(
        CP = z_contact,
        iCP = i_contact,
        delta = delta,
        noise = noise
      )
      return(append.afmdata(afmdata,CP))
    }
    }else{
      stop("No afmdata class input provided.")
    }
  }