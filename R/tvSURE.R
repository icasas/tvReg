#' Time-Varying Seemingly Unrelated Regression Equations Model
#'
#' Fits a set of balanced linear structural equations using Time-varying Ordinary Least 
#' Squares (tvOLS), Time-varying Seemingly Unrelated Regression (tvGLS), when the error 
#' variance-covariance matrix is known, or Time-varying Feasible Seemingly Unrelated 
#' Regression (tvFGLS), when the error variance-covariance matrix is unknown.
#'
#'
#' This function wraps up the kernel smoothing "tvOLS" and "tvGLS" estimators. The former is used when
#' equations are considered independent while the later assumes that the error term is correlated
#' amongst equations. This relation is given in matrix "Sigma" which is used in the estimation. When
#' "Sigma" is known, the estimates are calculated via the "tvGLS", and via the "tvFGLS" when "Sigma"
#' is unknown and must be estimated.
#'
#' Bandwidth selection is of great importance in kernel smoothing methodologies and it is done
#' automatically by cross-validation. One important aspect in the current packages is that the
#' bandwidth is selected independently for each equation and then the average is taken to use the
#' same bandwidth for each equation. It has been shown in Casas et al. (2017) that using
#' different bandwidths for each equation is in general a bad practice, even for uncorrelated equations.
#' Even though, the user may be able to use different bandwidths calling functions \code{\link{bw}} and
#' \code{\link{tvGLS}} separatedly.
#'
#' A system consists of "neq" number of equations with "obs" number of observations each and a number of
#' variables not necessarily equal for all equations. The matrix notation is:
#' \deqn{Y_{t} = X_t \beta_{t}+u_{t}}
#' where \eqn{Y_t  = (y_{1t}, y_{2t}, \ldots, y_{neq t})'}, \eqn{X_t = diag (x_{1t}, x_{2t}, \ldots, x_{neq t})}
#' and \eqn{\beta_{t} = \left(\beta _{1t}', \ldots, \beta _{neq t}'\right)'} is a vector of order the
#' total number of variables in the system. The error vector \eqn{u_{t} = (u_{1t}, u_{2t}, \ldots, u_{neq t})'}
#' has zero mean and  covariance matrix \eqn{E(u_t u_t') = \Sigma_t}.
#'
#' @seealso \code{\link{systemfit}}
#' @references
#' Casas, I., Ferreira, E., and Orbe, S. (2017) Time-Varying Coefficient Estimation 
#' in SURE Models: Application to Portfolio Management. Available at SSRN: 
#' https://ssrn.com/abstract=3043137
#' 
#' Chen, X. B., Gao, J., Li, D., and Silvapulle, P (2017) Nonparametric Estimation and 
#' Forecasting for Time-Varying Coefficient Realized Volatility Models.
#' \emph{Journal of Business \& Economic Statistics}, pp.1-13
#'
#' Granger, C. W (2008) Non-Linear Models: Where Do We Go Next - Time Varying
#' Parameter Models? \emph{Studies in Nonlinear Dynamics \& Econometrics}, 12, pp. 1-11.
#'
#' Kristensen, D (2012) Non-parametric detection and estimation of structural change.
#' \emph{Econometrics Journal}, 15, pp. 420-461.
#'
#' Orbe, S., Ferreira, E., and Rodriguez-Poo, J (2004) On the estimation and testing of
#' time varying constraints in econometric models, \emph{Statistica Sinica}.
#'
#'
#' @keywords time varying regression models, nonparametric statistics
#' @aliases tvsure-class tvsure
#' @rdname tvSURE
#' @importFrom systemfit systemfit
#' @param formula A list of formulas, one for each equation.
#' @param z A vector containing the smoothing variable.
#' @param bw An opcional scalar or vector of length the number of equations. It represents the bandwidth in
#' the estimation of trend coefficients. If NULL, it is selected by cross validation.
#' @param data A matrix or data frame containing variables in the formula.
#' @param method A character, a matrix of dimensions neq x neq or an array of dimensions obs x neq x neq, where
#' \code{obs} is the number of observations and \code{neq} is the number of equations.
#' If method = \code{identity} or \code{tvOLS} (default) then the method used is a time-varying OLS.
#' If method is a matrix (constant over time) or an array, then the \code{tvGLS} is called.
#' If method = \code{tvFGLS}, then the covariance matrix is estimated nonparametrically and the
#' estimation of the system is done as a whole.
#' @param Sigma A matrix of dimensions neq x neq or an array of dimensions neq x neq x obs
#' (neq = number of equations, obs = number of observations). It represents
#' the covariance matrix of the error term. Only necessary for method \code{tvGLS}.
#' @param bw.cov An optional scalar. It represents the bandwidth in the "lc" nonparametric estimation of the
#' time-varying covariance matrix. If NULL, it is selected by cross validation.
#' @param est The nonparametric estimation method, one of "lc" (default) for linear constant or "ll" for local linear.
#' @param tkernel The type of kernel used in the coefficients estimation method, one of Epanesnikov ("Epa") or "Gaussian".
#' @param singular.ok	Logical. If FALSE, a singular model is an error.
#' @param R An optional nrest x nvar x neq (nrest =  number of restrictions, nvar = number of variables in each equation,
#' neq = number of equations).
#' @param r An optional vector of length the number of restrictions. By default it contains zeros.
#' @param control list of control parameters.  The default is constructed by
#' the function \code{\link{tvsure.control}}.  See the documentation of
#' \code{\link{tvsure.control}} for details.
#' @param ... Other parameters passed to specific methods.
#' @return \code{tvSURE} returns a list of the class \code{tvsure} containing the results of the whole system, results of the estimation
#' and confidence instervals if chosen.
#' The object of class \code{tvsure} have the following components:
#' \item{tvcoef}{An array of dimension obs x nvar x neq (obs = number of observations, nvar = number of variables
#' in each equation, neq = number of equations in the system) with the time-varying coefficients estimates.}
#' \item{Lower}{If \code{level} non equal zero, an array of dimension obs x nvar x neq containing the confidence 
#' interval lower band.}
#' \item{Upper}{If \code{level} non equal zero, an array of dimension obs x nvar x neq containing the confidence 
#' interval upper band.}
#' \item{Sigma}{An array of dimension obs x neq x neq with the estimates of the errors covariance matrix.}
#' \item{fitted}{The fitted values.}
#' \item{residuals}{Estimation residuals.}
#' \item{x}{A list with the regressors data.}
#' \item{y}{A matrix with the dependent variable data.}
#' \item{z}{A vector with the smoothing variable.}
#' \item{bw}{Bandwidth of mean estimation.}
#' \item{obs}{Integer specifying the number of observations in each equation (balanced sample).}
#' \item{neq}{Integer specifying the number of equations.}
#' \item{nvar}{Vector of integers specifying the number of variables in each equation.}
#' \item{method}{Estimation method.}
#' \item{est}{Nonparemtric estimation methodology.}
#' \item{tkernel}{Kernel type.}
#' \item{bw.cov}{Bandwidht of Sigma estimation.}
#' \item{level}{Confidence interval range.}
#' \item{runs}{Number of bootstrap replications.}
#' \item{tboot}{Type of bootstrap.}
#' \item{BOOT}{List with all bootstrap replications of \code{tvcoef}, if done.}
#' \item{R}{Restrictions matrix.}
#' \item{r}{Restrictions vector.}
#' \item{formula}{Initial formula.}
#' \item{call}{Matched call.}
#' @seealso \code{\link{CI}}, \code{\link{plot}}
#' @examples
#' data("Kmenta", package = "systemfit")
#' eqDemand <- consump ~ price + income
#' eqSupply <- consump ~ price + farmPrice + trend
#' system <- list(demand = eqDemand, supply = eqSupply)
#'
#' ## OLS estimation
#' ols.fit <- systemfit::systemfit(system, method = "SUR", data = Kmenta)
#'
#' ## tvOLS estimation with the local linear estimator
#' tvols.fit <- tvSURE(system, data = Kmenta,  est = "ll")
#'
#' ## FGLS estimation - SURE estimation
#' fgls1.fit <- systemfit::systemfit(system, data = Kmenta, method = "SUR")
#'
#' ## tvFGLS estimation - tvSURE estimation
#' tvfgls1.fit <- tvSURE(system, data = Kmenta, method = "tvFGLS")
#'
#' ## iterative FGLS estimation - SURE estimation
#' fgls2.fit <- systemfit::systemfit(system, data = Kmenta, method = "SUR", maxit = 100)
#'
#' ## iterative tvFGLS estimation - SUR estimation using the local linear
#' tvfgls2.fit <- tvSURE(system, data = Kmenta, method = "tvFGLS",
#' control = list(tol = 0.001, maxiter = 100))
#'
#' ## Estimation with 2 restrictions
#' Rrestr <- matrix(0, 2, 7)
#' Rrestr[1, 3] <-  1
#' Rrestr[1, 7] <- -1
#' Rrestr[2, 2] <- -1
#' Rrestr[2, 5] <-  1
#' qrestr <- c(0, 0.5)
#'
#' tvfgls.rest <- tvSURE(system, data = Kmenta, method = "tvFGLS",
#' R = Rrestr, r = qrestr, bw = tvfgls1.fit$bw, bw.cov = tvfgls1.fit$bw.cov)
#'
#'@seealso \code{\link{bw}} for bandwidth calculation, \code{\link{tvGLS}} for nonparametric
#'estimator and \code{\link{CI}} for confidence intervals.
#'@export
tvSURE <- function (formula, z = NULL, bw = NULL, data,  method = c("tvOLS", "tvFGLS", "tvGLS"),  
                    Sigma = NULL, est = c("lc", "ll"), tkernel = c("Epa", "Gaussian"),
                    bw.cov = NULL, singular.ok = TRUE, R = NULL, r = NULL,
                    control = tvsure.control(...), ...)
{
  is.data <- inherits(data, c("data.frame", "matrix"))
  if(!is.data)
    stop("\nParameter 'data' should be enter and it should be a matrix or a data.frame.\n")
  if(class(formula) != "list")
    stop("\n'formula' must be a list of formulas. \n")
  if(!all(lapply(formula, class) == "formula"))
    stop("\nList 'x' must contain only objects of class 'formula'")
  neq <- length(formula)
  if(neq < 2)
    stop("\nThe list 'formula' should contain at least two equations for multivariate analysis.\n")
  if(!is.null (Sigma))
    if(any(is.na(Sigma)))
      stop("\nNAs in Sigma.\n")
  method <- match.arg (method)
  if(is.null(Sigma) & method == "tvGLS")
  {
    method <- "tvOLS"
    warning("\nSigma is NULL, tvOLS will be performed.\n")
  }
  est <- match.arg(est)
  tkernel <- match.arg(tkernel)
  if(est %in% c("lc", "ll"))
    est <- "lc"
  if(tkernel %in% c("Epa","Gaussian"))
    tkernel <- "Epa"
  nvar <- numeric(neq)
  if(is.null(names(formula)))
  {
    eq.names <- paste("eq", c(1:neq), sep = "")
  }
  else
  {
    eq.names <- names(formula)
    if(sum(regexpr(" |_", eq.names) != -1) > 0)
      stop("\nEquation labels may not contain blanks (' ') or underscores ('_')")
  }
  results <- list()
  results$call <- match.call()
  callNoDots <- match.call(expand.dots = FALSE)
  mf <- callNoDots[c(1, match("data", names(callNoDots), 0))]
  mf$na.action <- as.name("na.pass")
  mf[[1]] <- as.name("model.frame")
  y  <- NULL
  x  <- list()
  y.names <- NULL
  for(i in 1:neq)
  {
    mf.eq <-  mf
    mf.eq$formula <- formula[[i]]
    eval.mf <-  eval(mf.eq)
    terms <- attr(eval.mf, "terms")
    y <- cbind(y, stats::model.extract(eval.mf, "response"))
    y.names <- c(y.names, formula[[i]][[2]])
    x[[i]] <- stats::model.matrix(terms, eval.mf)
    nvar[i] <- ncol(x[[i]])
    if(is.null(colnames(x[[i]])))
      colnames(x[[i]]) <- paste("X", i, 1:nvar[i], sep = "")
  }
  eq.names <- names(formula)
  if(is.null(eq.names))
    eq.names <- paste("Eq", 1:neq, sep = "")
  names(x) <- eq.names
  colnames(y) <- y.names
  obs <- nrow(y)
  if(!is.null(R))
  {
    R <- as.matrix(R)
    if(ncol(R) != sum(nvar))
      stop("\nWrong dimension of R, it should have as many columns as variables 
           in the whole system. \n")
    if (is.null(r))
      r <- rep(0, nrow(R))
    else if (length(r) == 1)
      r <- rep(r, nrow(R))
    else if (length(r) != nrow(R) & length(r) != 1)
      stop("\nWrong dimension of r, it should be as long as the number of 
           rows in R. \n")
  }
  if (method == "identity" | method == "tvOLS")
  {
    if (is.null(bw))
      bw <- bw(x = x, y = y, z = z, est = est, tkernel = tkernel, 
               singular.ok = singular.ok)
    else
    {
      if (any(bw < 5/obs))
        stop("\nAt least one of your bw bandwidths is smaller than 5/obs, 
             please increase! \n")
      else if (any(is.na(bw)))
        stop("\nThe bandwidth cannot be a no number.\n")
    }
    result <- tvGLS(x = x, y = y, z = z, bw = bw, R = R, r = r, est = est, 
                    tkernel = tkernel)
    Sigma <- array(rep(crossprod(result$residuals)/ (obs - neq), obs), dim = c(neq, neq, obs))
  }
  else if(method == "tvFGLS")
  {
    if (is.null(bw))
      bw <- bw(x = x, y = y, z = z, est = est, tkernel = tkernel, singular.ok = singular.ok)
    else
    {
      if (any(bw<5/obs))
        stop("\nAt least one of your bw bandwidths is smaller than 5/obs, please increase! \n")
      else if (any (is.na(bw)))
        stop("\nThe bandwidth cannot be a no number.\n")
    }
    result <- tvGLS(x = x, y = y, z = z, bw = bw, R = R, r = r, est = est, tkernel = tkernel)
    bw.cov <- bwCov(x = result$residuals, tkernel = tkernel)
    Sigma <- tvCov(x = result$residuals, bw = bw.cov, tkernel = tkernel)
    result <- tvGLS(x = x, y = y, z = z, bw = bw, Sigma = Sigma, R = R, r = r,
                    est = est, tkernel = tkernel)
    itertemp <- 0
    tol <- control$tol
    maxiter <- control$maxiter
    tolold <- sum(result$tvcoef^2)
    tolnew <- 0
    while((abs(tolold-tolnew)>tol) && (itertemp <= maxiter))
    {
      tolold <- tolnew
      Sigma <- tvCov(bw = bw.cov, x = result$residuals, tkernel = tkernel)
      temp <- tvGLS(x = x, y = y, z = z, bw = bw, Sigma = Sigma, R = R, r = r,
                    est = est, tkernel = tkernel)
      tolnew <- sqrt(sum((result$tvcoef - temp$tvcoef)^2)/sum(result$tvcoef^2))
      result <- temp
      itertemp <- itertemp + 1
    }
  }
  else if(method == "tvGLS")
  {
    if(is.matrix(Sigma))
    {
      if(ncol(Sigma) != neq | nrow(Sigma) != neq)
        stop("\nWrong dimensions of Sigma. \n.")
      Sigma2 <- array(0, dim = c(neq, neq, obs))
      for (t in 1:obs)
        Sigma2[, , t] <- Sigma
      Sigma <- Sigma2
    }
    else if (is.array(Sigma))
    {
      dimensions <- dim(Sigma)
      if(dimensions[3] != obs | dimensions[2] != neq | dimensions[1] != neq)
        stop("\nWrong dimensions of Sigma. \n.")
    }
    else
      stop("\nSigma must be a matrix of dimensions neq x neq or an array of dimensions
           neq x neq x obs. \n")

    if (is.null(bw))
      bw <- bw(x = x, y = y, z = z, Sigma = Sigma, est = est, tkernel = tkernel)
    else
    {
      if (any(bw < 5/obs))
        stop("\nAt least one of your bw bandwidths is smaller than 5/obs,
             please increase! \n")
      else if (any (is.na(bw)))
        stop("\nThe bandwidth cannot be a no number.\n")
    }
    result <- tvGLS(x = x, y = y, z = z, bw = mean(bw), Sigma = Sigma, R = R, r = r,
                    est = est, tkernel = tkernel)
  }
  tvcoef <- result$tvcoef
  resid <- result$residuals
  fitted <- result$fitted
  if(length(bw) == 1)
    names(bw) <- "bw.mean"
  else
    names(bw) <- paste("bw.", colnames(y), sep = "")
  colnames(resid) <- paste(colnames(y), "eq", 1:neq, sep = "")
  colnames(fitted) <- colnames(resid)
  var.names <- NULL
  for(i in 1:neq)
    var.names <- c(var.names, paste(colnames(x[[i]]), ".eq", i, sep = ""))
  colnames(tvcoef) <- var.names
  result <- list(tvcoef =  tvcoef, Lower = NULL, Upper = NULL, Sigma = Sigma,
                 fitted = fitted, residuals = resid, x = x, y = y, z = z,  bw = bw, 
                 obs = obs, neq = neq, nvar = nvar, method = method, est =  est,
                 tkernel = tkernel, bw.cov = bw.cov, 
                 level = 0, runs = 0, tboot = NULL, BOOT = NULL,
                 R = R, r = r, formula = formula, call = match.call())
  class(result) <- "tvsure"
  return(result)
}


