##' Time-Varying Coefficients Linear Models
##'
##' \code{tvLM} is used to fit a time-varying coefficients linear model
##'
##' Models for \code{tvLM} are specified symbolically using the same formula
##' format than function \code{lm}. A typical model has the form \emph{response} ~ \emph{terms}
##' where response is the (numeric) response vector and terms is a series of terms which
##' specifies a linear predictor for response. A terms specification of the form
##' first + second indicates all the terms in first together with all the terms
##' in second with duplicates removed. A specification of the form first:second indicates
##' the set of terms obtained by taking the interactions of all terms in first with all
##' terms in second. The specification first*second indicates the cross of first and second.
##' This is the same as first + second + first:second.
##'
##' A formula has an implied intercept term. To remove this use either
##' y ~ x - 1 or y ~ 0 + x.
##'
##' @rdname tvLM
##' @aliases tvlm-class tvlm
##' @keywords time varying linear regression models, nonparametric statistics
##' @param formula An object of class formula.
##' @param z A vector with the smoothing variable.
##' @param data An optional data frame or matrix.
##' @param bw An opcional scalar. It represents the bandwidth in
##' the estimation of trend coefficients. If NULL, it is selected by cross validation.
##' @param est The nonparametric estimation method, one of "lc" (default) for linear constant
##'  or "ll" for local linear.
##' @param tkernel The type of kernel used in the coefficients estimation method,
##' one of Epanesnikov ("Epa") or "Gaussian".
##' @param singular.ok	Logical. If FALSE, a singular model is an error.
##'
##' @return An object of class \code{tvlm}
##' The object of class \code{tvlm} have the following components:
##' \item{tvcoef}{A matrix of dimensions}
##' \item{fitted}{The fitted values.}
##' \item{residuals}{Estimation residuals.}
##' \item{x}{A matrix with the regressors data.}
##' \item{y}{A vector with the dependent variable data.}
##' \item{z}{A vector with the smoothing variable.}
##' \item{bw}{Bandwidth of mean estimation.}
##' \item{est}{Nonparametric estimation methodology.}
##' \item{tkernel}{Kernel used in estimation.}
##' \item{level}{Confidence interval range.}
##' \item{runs}{Number of bootstrap replications.}
##' \item{tboot}{Type of bootstrap.}
##' \item{BOOT}{List with all bootstrap replications of \code{tvcoef}, if done.}
##' \item{call}{Matched call.}
##'
##' @examples
##' ## Simulate a linear process with time-varying coefficient
##' ## as functions of scaled time.
##'
##' tau <- seq(1:1000)/1000
##' beta <- data.frame(beta1 = sin(2 * pi * tau), beta2 = 2 * tau)
##' X1 <- rnorm(1000)
##' X2 <- rchisq(1000, df = 4)
##' error <- rt(1000, df = 10)
##' y <- apply(cbind(X1, X2) * beta, 1, sum) + error
##' data <-data.frame(y = y, X1 = X1, X2 = X2)
##'
##' ## Estimate coefficients with lm and tvLM for comparison
##'
##' coef.lm <- stats::lm(y ~ 0 + X1 + X2, data = data)$coef
##' model.tvlm <- tvLM(y ~ 0 + X1 + X2, data = data)
##'
##' ## Plot the estimates of beta1
##' plot(tau,beta[, 1], type = "l", main = "", ylab = expression(beta[1]),
##' xlab = expression(tau), ylim = range(beta[,1], model.tvlm$tvcoef[, 1]))
##' abline(h = coef.lm[1], col = 2)
##' lines(tau, model.tvlm$tvcoef[, 1], col=4)
##' legend("topright", c(expression(beta[1]), "lm", "tvlm"),
##' col = c(1, 2, 4), bty = "n", lty = 1)
##'
##' ##Obtain the 90% confidence interval of the coefficients
##' model.tvlm <- CI (model.tvlm, level = 0.90, runs = 50)
##' plot(model.tvlm)
##'
##' ## Simulate a linear process with time-varying coefficient
##' ## as functions of a random variable z
##'
##' set.seed (42)
##' z <- stats::arima.sim (n = 1000, list(ar = c(0.3, -0.2),
##' ma = c(-0.1, 0.2)), sd = sqrt(0.2))
##' beta <- data.frame(beta1 = sin(2 * pi * z), beta2 = 2 * z)
##' y <- apply(cbind(X1, X2) * beta, 1, sum) + error
##' data<-data.frame(y = y, X1 = X1, X2 = X2, z = z)
##'
##' coef.lm <- stats::lm(y ~ 0 + X1 + X2, data = data)$coef
##' model.tvlm2 <- tvLM(y ~ 0 + X1 + X2, z = z, data = data,
##' bw = 0.5, est = "ll")
##' model.tvlm2 <- CI(model.tvlm2)
##'
##' ##plot the estimates of beta1
##' sort.index <- sort.int(z, index.return = TRUE)$ix
##' plot(z[sort.index], beta[sort.index, 1], type = "l", main = "",
##' ylab = expression(beta[1]), xlab = expression(z[t]),
##' ylim = range(beta[,1], model.tvlm2$tvcoef[, 1]))
##' abline(h = coef.lm[1], col = 2)
##' lines(z[sort.index], model.tvlm2$tvcoef[sort.index, 1], col = 4)
##' legend("topleft", c(expression(beta[1]), "lm", "tvlm"),
##' col = c(1, 2, 4), bty = "n", lty = 1)
##'
##' @references 
##' Cai, Z., Li, Q., Park, J. Y. (2009) Functional-coefficient models for nonstationary 
##' time series data, \emph{Journal of Econometrics}, Volume 148, pp. 101-113.
##' 
##' Robinson, P. (1989) Nonparametric estimation of time-varying parameters.  In
##' Hackl, P., editor, \emph{Statistical Analysis and Forecasting of Economic Structural
##' Change}. Springer, Berlin.


##' @seealso \code{\link{bw}} for bandwidth selection, \code{\link{tvOLS}} for the
##' estimation procedure and \code{\link{CI}} for confidence intervals.
##' @keywords time-varying coefficients regression, nonparametric
##' @export

tvLM<-function (formula, z = NULL, data, bw = NULL, tkernel = "Epa", est = "lc",
                singular.ok = TRUE)
{
  if(tkernel != "Epa" & tkernel != "Gaussian")
    tkernel <- "Epa"
  if(est != "lc" & est != "ll")
    est <- "lc"
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- stats::model.response(mf, "numeric")
  if (stats::is.empty.model(mt))
    stop ("No regressors in the model. \n")
  else
  {
    x <- stats::model.matrix(mt, mf)
    if(is.null(bw))
      bw <- bw(x = x, y = y, z = z, tkernel = tkernel, est = est, singular.ok = singular.ok)
    results <- tvOLS (x = x, y = y, z = z, bw = bw, tkernel = tkernel, est = est,
                      singular.ok = singular.ok)
  }
  nvar <- ncol(x)
  xnames<-colnames(x)
  if(is.null(xnames))
    xnames <- paste("X", 1:nvar, collate="&", sep="")
  tvcoef <-results$tvcoef
  fitted <- results$fitted
  resid <- results$resid
  colnames(tvcoef) <- xnames
  result <- list(tvcoef = tvcoef, Lower = NULL, Upper = NULL, fitted = fitted,
                 residuals=resid, x = x, y = y, z = z, bw = bw, est = est, tkernel = tkernel,
                 singular.ok = singular.ok, level = 0, runs = 0, tboot = NULL, 
                 BOOT = NULL, call = cl)
  class(result) <- "tvlm"
  return(result)
}
