#' Time-Varying Coefficients Linear Models
#'
#' \code{tvLM} is used to fit a time-varying coefficients linear model
#'
#' Models for \code{tvLM} are specified symbolically using the same formula
#' format than function \code{lm}. A typical model has the form \emph{response} ~ \emph{terms}
#' where response is the (numeric) response vector and terms is a series of terms which
#' specifies a linear predictor for response. A terms specification of the form
#' first + second indicates all the terms in first together with all the terms
#' in second with duplicates removed. A specification of the form first:second indicates
#' the set of terms obtained by taking the interactions of all terms in first with all
#' terms in second. The specification first*second indicates the cross of first and second.
#' This is the same as first + second + first:second.
#'
#' A formula has an implied intercept term. To remove this use either
#' y ~ x - 1 or y ~ 0 + x.
#'
#' @rdname tvLM
#' @aliases tvlm-class tvlm
#' @keywords time varying linear regression models, nonparametric statistics
#' @param formula An object of class formula.
#' @param z A vector with the smoothing variable.
#' @param data An optional data frame or matrix.
#' @param bw An opcional scalar. It represents the bandwidth in
#' the estimation of trend coefficients. If NULL, it is selected by cross validation.
#' @param est The nonparametric estimation method, one of "lc" (default) for linear constant
#'  or "ll" for local linear.
#' @param tkernel The type of kernel used in the coefficients estimation method,
#' one of Epanesnikov ("Epa") or "Gaussian".
#' @param singular.ok	Logical. If FALSE, a singular model is an error.
#'
#' @return An object of class \code{tvlm}
#' The object of class \code{tvlm} have the following components:
#' \item{tvcoef}{A matrix of dimensions}
#' \item{fitted}{The fitted values.}
#' \item{residuals}{Estimation residuals.}
#' \item{x}{A matrix with the regressors data.}
#' \item{y}{A vector with the dependent variable data.}
#' \item{z}{A vector with the smoothing variable.}
#' \item{bw}{Bandwidth of mean estimation.}
#' \item{est}{Nonparametric estimation methodology.}
#' \item{tkernel}{Kernel used in estimation.}
#' \item{level}{Confidence interval range.}
#' \item{runs}{Number of bootstrap replications.}
#' \item{tboot}{Type of bootstrap.}
#' \item{BOOT}{List with all bootstrap replications of \code{tvcoef}, if done.}
#' \item{call}{Matched call.}
#' 
#' @seealso \code{\link{bw}}, \code{\link{tvAR}}, \code{\link{confint}}, 
#' \code{\link{plot}}, \code{\link{print}} and \code{\link{summary}}
#' @examples
#' ## Simulate a linear process with time-varying coefficient
#' ## as functions of scaled time.
#'
#' tau <- seq(0, 1, length.out = 200)
#' beta <- data.frame(beta1 = sin(2 * pi * tau), beta2 = 2 * tau)
#' X1 <- rnorm(200)
#' X2 <- rchisq(200, df = 4)
#' error <- rt(200, df = 10)
#' y <- apply(cbind(X1, X2) * beta, 1, sum) + error
#' data <-data.frame(y = y, X1 = X1, X2 = X2)
#'
#' ## Estimate coefficients with lm and tvLM for comparison
#'
#' coef.lm <- stats::lm(y ~ 0 + X1 + X2, data = data)$coef
#' model.tvlm <- tvLM(y ~ 0 + X1 + X2, data = data, bw = 0.2)
#'
#' ## Estimate coefficients of different realized variance models
#' data("RV")
#
#' ##Bollerslev et al. (2016) SHARQ model
#' SHARQ <- lm (RVt ~ RVt_1_pos + RVt_1_neg + I(RVt_1_pos * RQt_1_sqrt) +
#' I(RVt_1_neg * RQt_1_sqrt) + RVt_1_5 + RVt_1_22, data = tail(RV, 2000))
#' 
#' #Casas et al. (2018) tvSHARQ model
#' tvSHARQ <- tvLM (RVt ~ RVt_1_pos + RVt_1_neg + RVt_1_5 + RVt_1_22, 
#' z = tail(RV$RQt_1_sqrt, 2000), data = tail(RV, 2000), bw = 0.002)
#' 
#' boxplot(data.frame(tvSHARQ = tvSHARQ$tvcoef[,2],
#' SHARQ = SHARQ$coef[2]+ SHARQ$coef[4] * tail(RV$RQt_1_sqrt, 2000)),
#' main = expression (RV[t-1]^{"+"}), outline = FALSE)
#' boxplot(data.frame(tvSHARQ = tvSHARQ$tvcoef[,3], 
#' SHARQ = SHARQ$coef[3]+ SHARQ$coef[5] * tail(RV$RQt_1_sqrt, 2000)),
#' main = expression (RV[t-1]^{"-"}), outline = FALSE)
#'
#' @references 
#'  
#' Bollerslev, T., Patton, A. J. and Quaedvlieg, R. (2016) Exploiting the 
#' errors: A simple approach for improved volatility forecasting. 
#' \emph{Journal of Econometrics}, 192, 1-18.
#' 
#' Casas, I., Mao, X. and Vega, H. (2018) Reexamining financial and economic 
#' predictability with new estimators of realized variance and variance 
#' risk premium. Url= http://pure.au.dk/portal/files/123066669/rp18_10.pdf
#' 
#' @keywords time-varying coefficients regression, nonparametric
#' @export

tvLM<-function (formula, z = NULL, data, bw = NULL, est = c("lc", "ll"), 
                tkernel = c("Epa", "Gaussian"), singular.ok = TRUE)
{
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  if(!(tkernel %in% c("Epa","Gaussian")))
    tkernel <- "Epa"
  if(!(est %in% c("lc", "ll")))
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
