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
#' @param formula An object of class formula.
#' @param z A vector with the smoothing variable.
#' @param ez (optional) A scalar or vector with the smoothing estimation values. If 
#' values are included then the vector \code{z} is used.
#' @param data An optional data frame or matrix.
#' @param bw An opcional scalar. It represents the bandwidth in
#' the estimation of trend coefficients. If NULL, it is selected by cross validation. 
#' @param cv.block A positive scalar with the size of the block in leave one block out cross-validation.
#' By default 'cv.block=0' meaning leave one out cross-validation.
#' @param est The nonparametric estimation method, one of "lc" (default) for linear constant
#'  or "ll" for local linear.
#' @param tkernel A character, either "Triweight" (default), "Epa" or "Gaussian" kernel function.
#' @param singular.ok	Logical. If FALSE, a singular model is an error.
#'
#' @return An object of class \code{tvlm}
#' The object of class \code{tvlm} have the following components:
#' \item{coefficients}{A matrix of dimensions}
#' \item{fitted}{The fitted values.}
#' \item{residuals}{Estimation residuals.}
#' \item{x}{A matrix with the regressors data.}
#' \item{y}{A vector with the dependent variable data.}
#' \item{z}{A vector with the smoothing variable.}
#' \item{ez}{A vector with the smoothing estimation variable.}
#' \item{bw}{Bandwidth of mean estimation.}
#' \item{est}{Nonparametric estimation methodology.}
#' \item{tkernel}{Kernel used in estimation.}
#' \item{level}{Confidence interval range.}
#' \item{runs}{Number of bootstrap replications.}
#' \item{tboot}{Type of bootstrap.}
#' \item{BOOT}{List with all bootstrap replications of \code{coefficients}, if done.}
#' 
#' @seealso \code{\link{bw}}, \code{\link{tvAR}}, \code{\link{confint}}, 
#' \code{\link{plot}}, \code{\link{print}} and \code{\link{summary}}
#' @examples
#' ## Simulate a linear process with time-varying coefficient
#' ## as functions of scaled time.
#' set.seed(42)
#' tau <- seq(1:200)/200
#' beta <- data.frame(beta1 = sin(2*pi*tau), beta2= 2*tau)
#' X1 <- rnorm(200)
#' X2 <- rchisq(200, df = 4)
#' error <- rt(200, df = 10)
#' y <- apply(cbind(X1, X2)*beta, 1, sum) + error
#' data <- data.frame(y = y, X1 = X1, X2 = X2)
#' ## Estimate coefficients with lm and tvLM for comparison
#'
#' coef.lm <- stats::lm(y ~ 0 + X1 + X2, data = data)$coef
#' tvlm.fit <- tvLM(y ~ 0 + X1 + X2, data = data, bw = 0.29)
#'
#' ## Estimate coefficients of different realized variance models
#' data("RV")
#' RV2 <- head(RV, 2000)
#' ##Bollerslev t al. (2016) HARQ model
#' HARQ <- with(RV2, lm(RV ~ RV_lag + I(RV_lag * RQ_lag_sqrt) + RV_week + RV_month))
#' 
#' #Casas et al. (2018) TVHARQ model
#' TVHARQ <- with(RV2, tvLM (RV ~ RV_lag + RV_week + RV_month, z = RQ_lag_sqrt, 
#'                          bw = 0.0061))
#' boxplot(data.frame(TVHARQ = TVHARQ$coefficients[,2] * RV2$RV_lag,
#'                    HARQ = (HARQ$coef[2] + HARQ$coef[3] * RV2$RQ_lag_sqrt)*RV2$RV_lag),
#'                    main = expression (RV[t-1]), outline = FALSE)
#'                  
#' @references 
#'  
#' Bollerslev, T., Patton, A. J. and Quaedvlieg, R. (2016) Exploiting the 
#' errors: A simple approach for improved volatility forecasting. 
#' \emph{Journal of Econometrics}, 192, 1-18.
#' 
#' Casas, I., Mao, X. and Veiga, H. (2018) Reexamining financial and economic 
#' predictability with new estimators of realized variance and variance 
#' risk premium. Url= http://pure.au.dk/portal/files/123066669/rp18_10.pdf
#' 
#' @export

tvLM<-function (formula, z = NULL, ez = NULL, data, bw = NULL, cv.block = 0, est = c("lc", "ll"), 
                tkernel = c("Triweight", "Epa", "Gaussian"), singular.ok = TRUE)
{
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
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
    {
      cat("Calculating regression bandwidth... ")
      bw <- bw(x = x, y = y, z = z, cv.block = cv.block, est = est, tkernel = tkernel, 
               singular.ok = singular.ok)
      cat("bw = ", bw, "\n")
    }
    results <- tvOLS (x = x, y = y, z = z, ez = ez, bw = bw, tkernel = tkernel, est = est,
                      singular.ok = singular.ok)
  }
  nvar <- ncol(x)
  xnames <- colnames(x)
  yname <- attr(mt, "variables")[[2]]
  if(is.null(xnames))
    xnames <- paste0("X", 1:nvar, collate="&")
  coefficients <-results$coefficients
  colnames(coefficients) <- xnames
  result <- list(coefficients = coefficients, Lower = NULL, Upper = NULL, fitted = results$fitted,
                 residuals = results$resid, x = x, y = y, z = z, ez = ez, bw = bw, 
                 cv.block = cv.block, obs = length(y), est = est, tkernel = tkernel,
                 singular.ok = singular.ok, level = 0, runs = 0, 
                 tboot = NULL, BOOT = NULL, call = match.call())
  class(result) <- "tvlm"
  return(result)
}

