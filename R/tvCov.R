#' Time-varying Variance-Covariance Estimation
#'
#' Estimation of a time-varying/funcional coefficients variance-covariance matrix using the local constant or the local linear kernel
#' smoothing methodologies.
#' 
#' @references Aslanidis, N. and Casas, I (2013) Nonparametric correlation models for portfolio
#' allocation. \emph{Journal of Banking and Finance}, 37, 2268-2283
#'
#' @param x A matrix.
#' @param z A vector with the variable over which coefficients are smooth over.
#' @param ez (optional) A scalar or vector with the smoothing values. If 
#' values are not included then the vector \code{z} is used instead.
#' @param bw (Optional) A scalar.
#' @param cv.block A positive scalar with the size of the block in leave-one block-out cross-validation.
#' By default 'cv.block=0' meaning leave-one-out cross-validation.
#' @param est A character, either "lc" or "ll" for local constant or local linear.
#' @param tkernel A character, either "Triweight, "Epa" or "Gaussian" kernel functions.
#'
#' @return An array of dimension neq x neq x obs.
#' 
#' @seealso \code{\link{bwCov}}, \code{\link{tvCor}}
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
tvCov <- function(x, z = NULL, ez = NULL, bw = NULL, cv.block = 0, est = c("lc", "ll"), 
                  tkernel = c("Triweight", "Epa", "Gaussian"))
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
    bw <- bwCov(x, z = z, cv.block = abs(cv.block), est, tkernel)
  if(length(bw) > 1)
    bw <- stats::median(bw)
  if(!is.null(z))
  {
    if(length(z) != obs)
      stop("\nDimensions of 'x' and 'z' are not compatible\n")
    grid <- z
  }
  else
    grid <- (1:obs)/obs
  if (is.null(ez))
    ez <- grid
  for (t in 1:obs)
  {
    tau0 <- grid - ez[t]
    kernel.bw <- .kernel(x = tau0, bw = bw, tkernel = tkernel)
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
#' @param z A vector.
#' @param est A character, either "lc" or "ll" for local constant or local linear.
#' @param tkernel A character, either "Triweight" (default), "Epa" or "Gaussian" kernel function.
#'
#' @return A scalar with the mean squared error.
#' @keywords internal
#'
.tvCov.cv <- function(bw, x, z = NULL, cv.block = 0, est = c("lc", "ll"), tkernel = c("Triweight", "Epa", "Gaussian"))
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
  if(!is.null(z))
  {
    if(length(z) != obs)
      stop("\nDimensions of 'x' and 'z' are not compatible\n")
    grid <- z
  }
  else
    grid <- (1:obs)/obs
  for (t in 1:obs)
  {
    tau0 <- grid - grid[t]
    kernel.bw <- .kernel(x = tau0, bw = bw, tkernel = tkernel)
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


#' Time-varying Correlation Estimation
#'
#' Estimation of a time-varying/functional coefficients correlation matrix using the local constant or the local linear kernel
#' smoothing methodologies.
#' 
#'
#' @param x A matrix.
#' @param z A vector with the variable over which coefficients are smooth over.
#' @param ez (optional) A scalar or vector with the smoothing values. If 
#' values are not included then the vector \code{z} is used instead.
#' @param bw (optional) A scalar.
#' @param cv.block A positive scalar with the size of the block in leave-one block-out cross-validation.
#' By default 'cv.block=0' meaning leave-one-out cross-validation.
#' @param est A character, either "lc" or "ll" for local constant or local linear.
#' @param tkernel A character, either "Triweight, "Epa" or "Gaussian" kernel functions.
#'
#' @return An array of dimension neq x neq x obs.
#' 
#' @seealso \code{\link{tvCov}}
#'
#' @examples
#' ##Generate two independent (uncorrelated series)
#' y <- cbind(rnorm(100, sd = 4), rnorm(100, sd = 1))
#'
#' ##Estimation variance-variance matrix. If the bandwidth is unknown, it can
#' ##calculated with function bwCov()
#' Rho.hat <-  tvCor(y, bw = 1.4)
#' 
#' ##The first time estimate
#' print(Rho.hat[,,1])
#' ##The mean over time of all estimates
#' print(apply(Rho.hat, 1:2, mean))

#' ##Generate two dependent variables
#' y <- MASS::mvrnorm(n = 100, mu = c(0,0), Sigma = cbind(c(1, -0.5), c(-0.5, 4)))
#' 
#' ##Estimation variance-variance matrix
#' Rho.hat <-  tvCor(y, bw = 3.2)
#' ##The first time estimate
#' print(Rho.hat[,,1])
#'
#' @export tvCor
#'
tvCor <- function(x, z = NULL, ez = NULL, bw = NULL, cv.block = 0, est = c("lc", "ll"), 
               tkernel = c("Triweight", "Epa", "Gaussian"))
{
  correlation <- tvCov (x, z, ez, bw, cv.block, est, tkernel)
  obs <- dim(correlation)[3]
  ncol <- dim(correlation)[2]
  for (t in 1:obs)
    correlation[,,t] <- cov2cor(correlation[,,t])
  return(correlation)
}

