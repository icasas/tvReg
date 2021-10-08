#' Time-Varying Coefficients Panel Data Models
#'
#' Fits a balanced panel data model using the Time-Varying Pooled Ordinary Least 
#' Squares, the Time-Varying Random Effects and the Time-Varying Fixed Effects models. 
#' 
#' This function wraps up the kernel smoothing time-varying coefficient pooled, random effects
#' and fixed effects estimators. 
#'
#' Bandwidth selection is of great importance in kernel smoothing methodologies and it is done
#' automatically by cross-validation. 
#'
#' A panel data model consists of "neq" elements in the cross-sectional dimention and
#'  "obs" number of time observations for each cross-section. All variables are
#' the same for each equation which have common coefficients. 
#' 
#' @references
#' Casas, I., Gao, J., Peng B., and Xie, S. (2019). Modelling Time-Varying Income Elasticities 
#' of Health Care Expenditure for the OECD. 
#' Available at SSRN: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3262326
#' 
#' Sun, Y., Carrol, R.J and Li, D. (2009). Semiparametric Estimation of Fixed-Effects Panel Data 
#' Varying Coefficient Models. \emph{Advances in Econometrics}, 25, pp. 101-129.
#' 
#' @aliases tvplm-class tvplm
#' @rdname tvPLM
#' @param formula An object of class formula.
#' @param z A vector containing the smoothing variable.
#' @param ez (optional) A scalar or vector with the smoothing estimation values. If 
#' values are included then the vector \code{z} is used.
#' @param data An optional data frame or matrix.
#' @param index	Indicates the individual and time indexes. 
#' @param bw An opcional scalar. It represents the bandwidth in
#' the estimation of trend coefficients. If NULL, it is selected by cross validation. 
#' @param bw.cov An optional scalar. It represents the bandwidth in the "lc" nonparametric estimation of the
#' time-varying covariance matrix. If NULL, it is selected by cross validation for method \code{"random"}.
#' @param cv.block A positive scalar with the size of the block in leave one block out cross-validation.
#' By default 'cv.block=0' meaning leave one out cross-validation.
#' @param method A character with the choice of panel model/estimation method:
#' If \code{method = "pooling"} (default) then the data is pooled estimated with time-varying OLS. 
#' No individual or time effects are estimated
#' If \code{method = "random"} then individual effects are considered random and independent
#' of the regressors.
#' If \code{method = "within"} then individual effects which might be correlated with 
#' the regressors are estimated.
#' @param est The nonparametric estimation method, one of "lc" (default) for linear constant
#' @param tkernel A character, either "Triweight" (default), "Epa" or "Gaussian" kernel function.
#' @param control list of control parameters.  The default is constructed by
#' the function \code{\link{tvreg.control}}.  See the documentation of
#' \code{\link{tvreg.control}} for details.
#' @param ... Other parameters passed to specific methods.
#' @return \code{tvPLM} returns a list of the class \code{tvplm} containing the results of model, results of the estimation
#' and confidence instervals if chosen.
#' The object of class \code{tvplm} have the following components:
#' \item{coefficients}{An array of dimension obs x nvar x neq (obs = number of observations, nvar = number of variables
#' in each equation, neq = number of equations in the system) with the time-varying coefficients estimates.}
#' \item{Lower}{If \code{level} non equal zero, an array of dimension obs x nvar x neq containing the confidence 
#' interval lower band.}
#' \item{Upper}{If \code{level} non equal zero, an array of dimension obs x nvar x neq containing the confidence 
#' interval upper band.}
#' \item{fitted}{The fitted values.}
#' \item{residuals}{Estimation residuals.}
#' \item{x}{A list with the regressors data.}
#' \item{y}{A matrix with the dependent variable data.}
#' \item{z}{A vector with the smoothing variable.}
#' \item{ez}{A vector with the smoothing estimation values.}
#' \item{alpha}{A vector with the individual fixed effects, if chosen.}
#' \item{bw}{Bandwidth of mean estimation.}
#' \item{totobs}{Integer specifying the total number of observations.}
#' \item{neq}{Integer specifying the number of cross-section observations.}
#' \item{obs}{Integer specifying the number of time observations per cross-section.}
#' \item{nvar}{Number of variables.}
#' \item{method}{Estimation method.}
#' \item{est}{Nonparemtric estimation methodology.}
#' \item{tkernel}{Kernel type.}
#' \item{level}{Confidence interval range.}
#' \item{runs}{Number of bootstrap replications.}
#' \item{tboot}{Type of bootstrap.}
#' \item{BOOT}{List with all bootstrap replications of \code{coefficients}, if done.}
#' \item{formula}{Initial formula.}
#' \item{call}{Matched call.}
#' 
#' @seealso \code{\link{bw}}, \code{\link{confint}}, \code{\link{plot}}, 
#' \code{\link{print}} and \code{\link{summary}}
#' 
#' @examples
#' data(OECD)
#' ##TVPOLS estimation of the model
#' tvpols <- tvPLM(lhe~lgdp+pop65+pop14+public, index = c("country", "year"),
#'  data = OECD, method ="pooling", bw = 0.3)
#' \dontrun{
#' tvfe <- tvPLM(lhe~lgdp+pop65+pop14+public, index = c("country", "year"),
#'  data = OECD, method ="within", bw = 0.8)
#' tvre <- tvPLM(lhe~lgdp+pop65+pop14+public, index = c("country", "year"),
#'  data = OECD, method ="random", bw = 0.3)
#' }
#' @export
tvPLM <- function (formula, z = NULL, ez = NULL, data, index = NULL, bw = NULL, bw.cov = NULL, cv.block = 0, 
                   method = c("pooling", "random", "within"),  est = c("lc", "ll"), 
                   tkernel = c("Triweight", "Epa", "Gaussian"), control = tvreg.control(...), ...)
{
  is.panel <- inherits(data, c("data.frame", "matrix", "pdata.frame"))
  if(!is.panel)
    stop("\nArgument 'data' should be entered and it should be a 'matrix', a 'data.frame' or 'pdata.frame'.\n")
  if (!inherits(formula, c("formula"))) 
      stop("\nArgument 'formula' should be entered of class 'formula'.")
  if (!inherits(data, "pdata.frame")) 
    data <- plm::pdata.frame(data, index)
  else
    if(dim(attr(data, "index"))[2] > 2)
      stop("'index' can be of length 2 at the most (one index variable for individual, time)")
  method <- match.arg (method)
  if(attr(terms(formula), "intercept"))
  {
    formula <- stats::update(formula,  ~ . + 0)
  }
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf$formula <- data
  mf$data <- formula
  data <- eval(mf, parent.frame())
  neq <- plm::pdim(data)$nT$n
  obs <- plm::pdim(data)$nT$T
  terms <- attr(data, "terms")
  y <- stats::model.extract(data, "response")
  if (stats::is.empty.model(data))
    stop ("\nNo regressors in the model. \n")
  if(length(y) != neq*obs)
    stop("\nNAs in the dependent variable. \n")
  x <- stats::model.matrix(terms, data)
  if(dim(x)[1] != neq*obs)
    stop("\nNAs in the regressors. \n")
  var.names <- colnames(x)
  nvar <- length(var.names)
  if(is.null(bw))
  {
    cat("Calculating regression bandwidth... ")
    bw <- bw(data, z = z, method = method, cv.block = cv.block, 
             est = est, tkernel = tkernel)
    cat("bw = ", bw, "\n")
  }
  if(method != "within")
  {
    result <- tvRE (x = x, y = y, z = z, ez = ez, bw = bw, Sigma = NULL,
                      neq = neq, obs = obs,  est = est, tkernel = tkernel)
    if(method == "random")
    {
      result <- tvRE(x = x, y = y, z = z, ez = ez, bw = bw, Sigma = result$Sigma, 
                     neq = neq, obs = obs, est = est, tkernel = tkernel)
      itertemp <- 1
      tol <- control$tol
      maxiter <- control$maxiter
      tolold <- sum(result$coefficients^2)
      tolnew <- 0
      while((abs(tolold-tolnew)>tol) && (itertemp < maxiter))
      {
        tolold <- tolnew
        temp <- tvRE(x = x, y = y, z = z, ez = ez, bw = bw, Sigma = result$Sigma, 
                               neq = neq, obs = obs, est = est, tkernel = tkernel)
        tolnew <- sqrt(sum((result$coefficients - temp$coefficients)^2)/sum(result$coefficients^2))
        result <- temp
        itertemp <- itertemp + 1
      }
    }
  }
  else 
    result <- tvFE (x = x, y = y, z = z, ez = ez, bw = bw,
                    neq = neq, obs = obs, est = est, tkernel = tkernel)
  alpha <- result$alpha
  coefficients <- result$coefficients
  resid <- result$residuals
  fitted <- result$fitted
  colnames(coefficients) <- var.names
  result <- list(coefficients =  coefficients, alpha = alpha, Lower = NULL, Upper = NULL, 
                 fitted = fitted, residuals = resid, x = x, y = y, z = z, ez = ez, bw = bw,
                 Sigma = result$Sigma, cv.block = cv.block, neq = neq, obs = obs, nvar = nvar, 
                 index = plm::index(data), method = method, est =  est, tkernel = tkernel, 
                 control = control, level = 0, runs = 0, tboot = NULL, BOOT = NULL,
                 formula = formula, call = match.call())
  class(result) <- "tvplm"
  return(result)
}


