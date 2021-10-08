#' Time-varying Variance-Covariance Estimation
#'
#' Estimation of a time-varying variance-covariance matrix using the local constant or the local linear kernel
#' smoothing methodologies.
#' 
#' @references Aslanidis, N. and Casas, I (2013) Nonparametric correlation models for portfolio
#' allocation. \emph{Journal of Banking \& Finance}, 37, 2268-2283
#'
#' @param x A matrix.
#' @param bw A scalar.
#' @param cv.block A positive scalar with the size of the block in leave-one block-out cross-validation.
#' By default 'cv.block=0' meaning leave-one-out cross-validation.
#' @param est A character, either "lc" or "ll" for local constant or local linear.
#' @param tkernel A character, either "Triweight, "Epa" or "Gaussian" kernel functions.
#'
#' @return A matrix of dimension obs x neq x neq.
#' 
#' @seealso \code{\link{bwCov}}
#'
#' @examples
#' ##Generate two independent (uncorrelated series)
#' y <- cbind(rnorm(100, sd = 4), rnorm(100, sd = 1))
#'
#' ##Estimation variance-variance matrix. If the bandwidth is unknown, it can
#' ##calculated with function bwCov()
#' Sigma.hat <-  tvCov(y, bw = 1.4)
#' 
#' ##The first time estimate
#' print(Sigma.hat[,,1])
#' ##The mean over time of all estimates
#' print(apply(Sigma.hat, 1:2, mean))

#' ##Generate two dependent variables
#' y <- MASS::mvrnorm(n = 100, mu = c(0,0), Sigma = cbind(c(1, -0.5), c(-0.5, 4)))
#' 
#' ##Estimation variance-variance matrix
#' Sigma.hat <-  tvCov(y, bw = 3.2)
#' ##The first time estimate
#' print(Sigma.hat[,,1])
#'
#' @export tvCov
#'
tvCov <- function(x, bw = NULL, cv.block = 0, est = c("lc", "ll"), tkernel = c("Triweight", "Epa", "Gaussian"))
{
  if(!inherits(x, c("matrix", "data.frame")))
    stop("'x' should be a matrix or a data.frame.\n")
  x <- as.matrix(x)
  obs <- NROW(x)
  neq <- NCOL(x)
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  Sigma <- array(0, dim = c(neq, neq, obs))
  resid.2 <- numeric(obs)
  if(is.null(bw))
    bw <- bwCov(x, cv.block = abs(cv.block), est, tkernel)
  if(length(bw) > 1)
    bw <- stats::median(bw)
  time.grid <- 1:obs
  for (t in 1:obs)
  {
    tau0 <- (time.grid - t)/obs
    kernel.bw <- .kernel(tau0, bw, tkernel = tkernel)/bw
    w0 <- sum(kernel.bw)
    r.2 <- matrix(0, neq, neq)
    if(est == "lc")
    {
      for (index in which(kernel.bw != 0))
      {
        r.2 <- tcrossprod(x[index,]) * kernel.bw[index]/w0 + r.2
      }
    }
    else if (est == "ll")
    {
      w1 <- sum(kernel.bw * tau0)
      w2 <- sum(kernel.bw * (tau0)^2)
      den.c <- w0 * w2 - w1^2
      if (den.c == 0) return(.Machine$double.xmax)
      for (index in which(kernel.bw != 0))
      {
        r.2 <- tcrossprod(x[index,]) * kernel.bw[index] * (w2-w1 * tau0[index])/den.c + r.2
      }
    }
    Sigma[,,t] <- r.2
  }
  return(Sigma)
}

#' @name tvReg-internals
#' @param x A matrix.
#' @param bw A scalar.
#' @param est A character, either "lc" or "ll" for local constant or local linear.
#' @param tkernel A character, either "Triweight" (default), "Epa" or "Gaussian" kernel function.
#'
#' @return A scalar with the mean squared error.
#' @keywords internal
#'
.tvCov.cv <- function(bw, x, cv.block = 0, est = c("lc", "ll"), tkernel = c("Triweight", "Epa", "Gaussian"))
{
  x <- as.matrix(x)
  obs <- NROW(x)
  neq <- NCOL(x)
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  Sigma <-  array(0, dim = c(neq, neq, obs))
  resid.2 <- numeric(obs)
  if(length(bw)>1) 
    bw <- stats::median(bw)
  grid <- 1:obs/obs
  for (t in 1:obs)
  {
    tau0 <- grid - grid[t]
    kernel.bw <- .kernel(tau0, bw, tkernel = tkernel)/bw
    kernel.bw[max(1, t-cv.block):min(t+cv.block, obs)] <- 0
    w0 <- sum(kernel.bw)
    k.index <- which(kernel.bw != 0)
    if (sum(k.index != 0) < 3)
      return (.Machine$double.xmax)
    r.2 <- matrix(0, neq, neq)
    if(est == "lc")
    {
      for (index in k.index)
      {
        r.2 <- tcrossprod(x[index,]) * kernel.bw[index]/w0 + r.2
      }
    }
    else if (est == "ll")
    {
      w1 <- sum(kernel.bw * tau0)
      w2 <- sum(kernel.bw * (tau0)^2)
      den.c <- w0 * w2 - w1^2
      if (den.c == 0) return(.Machine$double.xmax)
      for (index in k.index)
      {
        r.2 <- tcrossprod(x[index,]) * kernel.bw[index] * (w2 - w1 * tau0[index])/den.c + r.2
      }
    }
    Sigma[,,t] <- r.2
    resid.2[t] <- sum((tcrossprod(x[t,]) - Sigma[,,t])^2)
  }
  return(mean(resid.2))
}

#' Estimated variance-covariance matrix from TVPOLS
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
