#' Time-Varying Autoregressive Model
#'
#' \code{tvAR} is used to fit an autorregressive model with time varying coefficients.
#' 
#' It is a special case of linear model in which the regressors are lags of the
#' dependent variable. If any variable is included in the \code{xreg} term, these are added 
#' to the regressors matrix. A time-varying coefficients linear regression (with an
#' intercept if type = "const") is fitted.
#'
#' @references 
#' 
#' Cai, Z. (2007) Trending time-varying coefficient time series with serially
#' correlated errors, \emph{Journal of Econometrics}, 136, pp. 163-188.
#' 
#' @param y A vector with the dependent variable.
#' @param p A scalar indicating the number of lags in the model.
#' @param bw An opcional scalar or vector of length the number of equations. It represents
#' the bandwidth in the estimation of coefficients. If NULL, it is selected
#' by cross validation.
#' @inheritParams tvLM
#' @param fixed (optional) numeric vector of the same length as the total number of parameters.
#' If supplied, only NA entries in fixed will be varied.
#' 
#' @return An object of class \code{tvar} with the following components:
#' \item{tvcoef}{A vector of dimension obs (obs = number of observations - number lags),
#'  with the time-varying coefficients estimates.}
#' \item{fitted}{The fitted values.}
#' \item{residuals}{Estimation residuals.}
#' \item{x}{A matrix of model data, with lagged y and exogenous variables.}
#' \item{y}{A vector with the dependent data used in the model.}
#' \item{z}{A vector with the smoothing variable in the model.}
#' \item{ez}{A vector with the smoothing estimation values.}
#' \item{y.orig}{A vector with the original variable y.}
#' \item{bw}{Bandwidth of mean estimation.}
#' \item{type}{Whether the model has a constant or not.}
#' \item{exogen}{A matrix or data.frame with other exogenous variables.}
#' \item{p}{Number of lags}
#' \item{obs}{Number of observations in estimation.}
#' \item{totobs}{Number of observations in the original set.}
#' \item{level}{Confidence interval range.}
#' \item{runs}{Number of bootstrap replications.}
#' \item{tboot}{Type of bootstrap.}
#' \item{BOOT}{List with all bootstrap replications of \code{tvcoef}, if done.}
#' 
#' @seealso  \code{\link{bw}}, \code{\link{tvLM}}, \code{\link{confint}}, 
#' \code{\link{plot}}, \code{\link{print}} and \code{\link{summary}}
#' 
#' @examples
#'
#' ## Estimate coefficients of different realized variance models
#' data("RV")
#' RV2 <- head(RV, 2000)
#' RV <- RV2$RV
#' RV_week <- RV2$RV_week
#' RV_month <- RV2$RV_month
#' RQ <- RV2$RQ_lag_sqrt
#
#' ##Corsi (2009) HAR model
#' HAR <- arima(RV, order = c(1, 0, 0), xreg = cbind (RV_week, RV_month))
#' print(HAR)
#' 
#' ##Chen et al (2017) TVC-HAR model 
#' TVCHAR <- tvAR (RV, p = 1, exogen = cbind (RV_week, RV_month), bw = 20)
#' print(TVCHAR)
#' 
#' ##Casas et al (2018) TV-HARQ model
#' tvHARQ <- tvAR (RV, p = 1, exogen = cbind (RV_week, RV_month), 
#' z=RQ, bw = 0.0062)
#' print(tvHARQ)
#' 
#' @aliases tvar-class tvar
#' @references
#' 
#' Casas, I., Mao, X. and Veiga, H. (2018) Reexamining financial and economic 
#' predictability with new estimators of realized variance and variance 
#' risk premium. Url= http://pure.au.dk/portal/files/123066669/rp18_10.pdf
#'
#' Chen, X. B., Gao, J., Li, D., and Silvapulle, P (2017) Nonparametric estimation and 
#' forecasting for time-varying coefficient realized volatility models,
#' \emph{Journal of Business \& Economic Statistics}, online, 1-13.
#' 
#' Corsi, F. (2009) A simple approximate long-memory model of realized 
#' volatility. \emph{Journal of Financial Econometrics}, 7, 174-196.
#'
#' @rdname tvAR
#' @inheritParams tvVAR
#' @export
tvAR <- function (y, p = 1, z = NULL, ez = NULL, bw = NULL, cv.block = 0, type = c("const", "none"), exogen = NULL,
                  fixed = NULL, est = c("lc", "ll"), tkernel = c("Epa", "Gaussian"), singular.ok = TRUE)
{
  if (any(is.na(y)))
    stop("\nNAs in y.\n")
  if (p < 1)
    stop("p should be a positive number. \n")
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  type <- match.arg(type)
  y.orig <- y
  obs <- length(y)
  sample <- obs - p
  ylags <- as.matrix(stats::embed(y, dimension = p + 1))
  colnames(ylags) <- c("y", paste("y.l", 1:p, sep=""))
  yend <- ylags[, 1]
  ylags <- ylags[, -1, drop = FALSE]
  rhs <- ylags
  colnames(rhs) <- make.names(colnames(ylags))
  if (type == "const") {
    rhs <- cbind(rep(1, sample), ylags)
    colnames(rhs) <- c("(Intercept)", colnames(ylags))
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
  colnames(datamat) <- colnames(rhs)[mask]
  if(is.null(bw))
  {
    cat("Calculating regression bandwidth...\n")
    bw <- bw(x = datamat, y = yend, z = z, tkernel = tkernel, est = est, singular.ok = singular.ok)
  }
  results <- tvOLS(x = datamat, y = yend, z = z, ez = ez, bw = bw, est = est, tkernel = tkernel,
                   singular.ok = singular.ok)
  tvcoef <- results$tvcoef
  colnames(tvcoef) <- colnames(rhs)[mask]
  result <- list(tvcoef = tvcoef, Lower = NULL, Upper = NULL,  fitted = results$fitted,
                 residuals = results$resid, x = datamat, y = yend, z = z, ez = ez, y.orig = y.orig, 
                 bw = bw, cv.block = cv.block, mask = mask, exogen = exogen, p = p, type = type, obs = sample, 
                 totobs = sample + p, est = est, tkernel = tkernel, level = 0, runs = 0, tboot = NULL, BOOT = NULL)
  class(result) <- c("tvar")
  return(result)
}


