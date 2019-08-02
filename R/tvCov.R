#' Time-varying Variance-Covariance Estimation
#'
#' Estimation of a time-varying variance-covariance matrix using the local constant or the local linear kernel
#' smoothing methodologies.
#' 
#' @references Aslanidis, N. and Casas, I (2013) Nonparametric correlation models for portfolio
#' allocation. \emph{Journal of Banking \& Finance}, 37, 2268-2283
#'
#' @importFrom MASS mvrnorm
#' @param x A matrix.
#' @param bw A scalar.
#' @param est A character, either "lc" or "ll" for local constant or local linear.
#' @param tkernel A character, either "Gaussian" or "Epa" kernel types.
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
tvCov <- function(x, bw, est = c("lc", "ll"), tkernel = c("Epa", "Gaussian"))
{
  x <- as.matrix(x)
  obs <- NROW(x)
  neq <- NCOL(x)
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  Sigma <- array(0, dim = c(neq, neq, obs))
  resid.2 <- numeric(obs)
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
#' @param tkernel A character, either "Gaussian" or "Epa" kernel types.
#'
#' @return A scalar with the mean squared error.
#' @keywords internal
#'
.tvCov.cv <- function(bw, x, cv.block = 0, est = c("lc", "ll"), tkernel = c("Epa", "Gaussian"))
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
  time.grid <- 1:obs
  for (t in 1:obs)
  {
    tau0 <- (time.grid-t)/obs
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
