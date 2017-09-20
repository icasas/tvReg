##' Time-Varying Autoregressive Model
##'
##' \code{tvAR} is used to fit an autorregressive model with time varying coefficients.
##' 
##' It is a special case of linear model in which the regressors are lags of the
##' dependent variable. If any variable is included in the \code{xreg} term, these are added 
##' to the regressors matrix. A time-varying coefficients linear regression (with an
##' intercept if type = "const") is fitted.
##'
##' @references 
##' 
##' Cai, Z. (2007) Trending time-varying coefficient time series with serially
##' correlated errors, \emph{Journal of Econometrics}, 136, pp. 163-188.
##' 
##' @param y A vector with the dependent variable.
##' @param p A scalar indicating the number of lags in the model.
##' @param z A vector with the smoothing variable.
##' @param bw An opcional scalar or vector of length the number of equations. It represents
##' the bandwidth in the estimation of coefficients. If NULL, it is selected
##' by cross validation.
##' @inheritParams tvLM
##' @param fixed (optional) numeric vector of the same length as the total number of parameters.
##' If supplied, only NA entries in fixed will be varied.
##'
##' @keywords time-varying coefficients, nonparametric statistics
##' @aliases tvar-class tvar
##'
##' @return An object of class 'tvar'
##' The object of class \code{tvar} have the following components:
##' \item{tvcoef}{A vector of dimension obs (obs = number of observations - number lags),
##'  with the time-varying coefficients estimates.}
##' \item{fitted}{The fitted values.}
##' \item{residuals}{Estimation residuals.}
##' \item{x}{A matrix of model data, with lagged y and exogenous variables.}
##' \item{y}{A vector with the dependent data used in the model.}
##' \item{z}{A vector with the smoothing variable in the model.}
##' \item{bw}{Bandwidth of mean estimation.}
##' \item{type}{Whether the model has a constant or not.}
##' \item{exogen}{A matrix or data.frame with other exogenous variables.}
##' \item{p}{Number of lags}
##' \item{obs}{Number of observations in estimation.}
##' \item{totobs}{Number of observations in the original set.}
##' \item{level}{Confidence interval range.}
##' \item{runs}{Number of bootstrap replications.}
##' \item{tboot}{Type of bootstrap.}
##' \item{BOOT}{List with all bootstrap replications of \code{tvcoef}, if done.}
##' \item{call}{Matched call.}
##'
##' @examples
##' ## Simulate an tvAR(2) process
##'
##' tt <- (1:1000)/1000
##' beta <- cbind( 0.5 * cos (2 * pi * tt), (tt - 0.5)^2)
##' y <- numeric(1000)
##' y[1] <- 0.5
##' y[2] <- -0.2
##'
##' ## y(t) = beta1(t) y(t-1) + beta2(t) y(t-2) + ut
##'
##' for (t in 3:1000)
##' {
##'   y[t] <- y[(t-1):(t-2)] %*% beta[t,] + rnorm(1)
##' }
##' Y <- tail (y, 500)
##'
##' ## Estimate coefficients of process Y with ar.ols and tvAR
##' ## and compare them in a plot
##'
##' tvAR.2p <- tvAR(Y, p = 2, type = "none", est = "ll")
##' AR.2p <- ar.ols(Y, aic = FALSE, order = 2, intercept = FALSE, demean = FALSE )
##' plot(tail(beta[, 1], 500), ylim=range(tvAR.2p$tvcoef[, 1], tail(beta[, 1], 500)),
##' xlab = "", ylab = "", cex = 0.5, pch = 20)
##' abline(h = AR.2p$ar[1], col = 2)
##' lines(tvAR.2p$tvcoef[, 1], col = 4)
##' legend("topleft", c(expression(beta[1]),"AR", "tvAR"), col = c(1, 2, 4),
##' lty = 1, bty = "n")
##'
##' ## Estimate only coefficient from odd lags and the intercept
##' tvAR.6p <- tvAR(Y, p = 6, type = "const",
##' fixed = c(NA, 0, NA, 0, NA, 0, NA), est = "ll")
##' 
##' ##' ## Generation of model with coefficients depending of a random variable
##' z <- arima.sim(n = 1000, list(ma = c(-0.2279, 0.2488)))
##' beta <- (z - 0.5)^2
##' y <- numeric(1000)
##' y[1] <- -1
##' 
##' ##y(t) = beta1(z(t)) y(t-1) + ut
##'  
##' for (t in 3:1000)
##' {
##'   y[t] <- y[(t-1)] %*% beta[t] + rnorm(1)
##' }
##' Y <- tail (y, 500)
##' Z <- tail(z, 500)
##' 
##' ## Estimate coefficients of process Y with ar.ols and tvAR
##' ## and compare them in a plot
##' 
##' tvAR.2p.z <- tvAR(Y, z = Z, p = 1, type = "none", est = "ll")
##' AR.2p <- ar.ols(Y, aic = FALSE, order = 1, intercept = FALSE, demean = FALSE )
##' index <- sort.int (tvAR.2p.z$z, index.return = TRUE)$ix
##' z.sort <- tvAR.2p.z$z[index]
##' beta <- tail(beta, length(z.sort))[index]
##' beta.hat <- tvAR.2p.z$tvcoef[index]
##' plot(z.sort, beta, ylim=range(beta.hat, beta),
##'    xlab = "", ylab = "", type="l")
##' abline(h = AR.2p$ar[1], col = 2)
##' lines(z.sort, beta.hat, col = 4)
##' legend("top", c(expression(beta),"AR", "tvAR"), col = c(1, 2, 4),
##'        lty = 1, bty = "n")
##'
##' @seealso \code{\link{tvLM}} for estimation of time-varying coefficients linear models,
##'  and \code{\link{CI}} for confidence intervals.
##' @rdname tvAR
##' @inheritParams tvVAR
##' @export
tvAR <- function (y, p = 1, z = NULL, bw = NULL, type = c("const", "none"), exogen = NULL,
                  fixed = NULL, tkernel = "Epa", est = "lc", singular.ok = TRUE)
{
  if (any(is.na(y)))
    stop("\nNAs in y.\n")
  if (p < 1)
    stop("p should be a positive number. \n")
  if(tkernel != "Epa" & tkernel != "Gaussian")
    tkernel <- "Epa"
  if(est != "lc" & est != "ll")
    est <- "lc"
  y.orig <- y
  type <- match.arg(type)
  obs <- length(y)
  sample <- obs - p
  ylags <- as.matrix(stats::embed(y, dimension = p + 1))
  colnames(ylags) <- c("y", paste("y.l", 1:p, sep=""))
  yend <- ylags[, 1]
  ylags <- ylags[, -1, drop = FALSE]
  rhs <- ylags
  colnames(rhs) <- make.names(colnames(ylags))
  if (type == "const") {
    rhs <- cbind( ylags, rep(1, sample))
    colnames(rhs) <- c(colnames(ylags), "Intercept")
  }
  if (!(is.null(exogen))) 
  {
    exogen <- as.matrix(exogen)
    if (!identical(NROW(exogen), NROW(y)))
      stop("\nDifferent row size of 'y' and exogen.\n")
    if (is.null(colnames(exogen))) 
      colnames(exogen) <- paste("exo", 1:ncol(exogen),
                                sep = "")
    colnames(exogen) <- make.names(colnames(exogen))
    tmp <- colnames(rhs)
    rhs <- cbind(rhs, exogen[-c(1:p), ])
    colnames(rhs) <- c(tmp, colnames(exogen))
  }
  if (!is.null(z))
  {
    if(!identical(NROW(z), NROW(y)))
      stop("\nDifferent row size of 'y' and 'z'.\n")
    z <- z[-c(1:p)]
  }
  nvar <- ncol(rhs)
  if (is.null(fixed))
    fixed <- rep(NA_real_, nvar)
  else if (length(fixed) != nvar)
    stop("\nWrong length for 'fixed'\n")
  mask <- is.na(fixed)
  datamat <- as.matrix(rhs[, mask])
  if(is.null(bw))
    bw <- bw(x = datamat, y = yend, z = z, tkernel = tkernel, est = est, singular.ok = singular.ok)
  results <- tvOLS(x = datamat, y = yend, z = z, bw = bw, est = est, tkernel = tkernel,
                   singular.ok = singular.ok)
  tvcoef <- results$tvcoef
  colnames(tvcoef) <- colnames(rhs)[mask]
  result <- list(tvcoef = tvcoef, Lower = NULL, Upper = NULL,  fitted = results$fitted,
                 residuals = results$resid, x = datamat, y = yend, z = z, y.orig = y.orig, 
                 mask = mask, exogen = exogen, p = p, type = type, obs = sample, 
                 totobs = sample + p, est = est, tkernel = tkernel, bw = bw, level = 0,
                 runs = 0, tboot = NULL, BOOT = NULL, call = match.call())
  class(result) <- c("tvar", "tvlm")
  return(result)
}


