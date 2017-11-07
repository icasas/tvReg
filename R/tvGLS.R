#' Time-varying Generalised Least Squares
#'
#' \code{tvGLS} estimates time-varying coefficients of SURE using the kernel smoothing GLS.
#'
#'
#' @param x an object used to select a method.
#' @param ... Other parameters passed to specific methods.
#' @return \code{tvGLS} returns a list containing:
#' \item{tvcoef}{An array of dimension obs x nvar x neq (obs = number of observations, nvar = number of variables
#' in each equation, neq = number of equations in the system) with the time-varying coefficients estimates.}
#' \item{fited}{A matrix of dimension obs x neq with the fited values from the estimation.}
#' \item{residuals}{A matrix of dimension obs x neq with the residuals from the estimation.}
#' @export
#' @import Matrix

tvGLS<- function(x, ...) UseMethod("tvGLS", x)

#' Estimate Time-varying Coefficients
#'
#' \code{tvGLS} is used to estimate time-varying coefficients SURE using the kernel smoothing
#' generalised least square.
#'
#' The classical GLS estimator must be modified to generate a set of coefficients changing over time.
#' The \code{tvGLS} finds a GLS estimate at a given point in time \emph{t} using the data near by.
#' The size of the data window used is given by the bandwidth. The closest a point is to \emph{t},
#' the larger is its effect on the estimation which is given by the kernel. In this programme,
#' the two possible kernels are the Epanechnikov and Gaussian. As in the classical GLS, the covariance
#' matrix is involved in the estimation formula. If this matrix is NULL or the identity, then the
#' programme returns the OLS estimates for time-varying coefficients.
#'
#' Note, that unless with the tvSURE, the tvGLS may run with one common bandwidth for all
#' equations or with a different bandwidths for each equation.
#'
#' @rdname tvGLS
#' @param y A matrix.
#' @param z A vector with the variable over which coefficients are smooth over.
#' @param bw A numeric vector.
#' @param Sigma An array.
#' @param R A matrix.
#' @param r A numeric vector.
#' @param est Either "lc" or "ll".
#' @param tkernel Either "Gaussian" or "Epa".
#' @return A list with the estimates, fitted and residuals values.
#'
#' @examples
#' data(FF5F)
#' x <- list()
#' ## SMALL/LoBM porfolios time-varying three factor model
#' x[[1]] <- FF5F[, c("NA.Mkt.RF", "NA.SMB",  "NA.HML", "NA.RMW", "NA.CMA")]
#' x[[2]] <- FF5F[, c("JP.Mkt.RF", "JP.SMB",  "JP.HML", "JP.RMW", "JP.CMA")]
#' x[[3]] <- FF5F[, c("AP.Mkt.RF", "AP.SMB",  "AP.HML", "AP.RMW", "AP.CMA")]
#' x[[4]] <- FF5F[, c("EU.Mkt.RF", "EU.SMB",  "EU.HML", "EU.RMW", "EU.CMA")]
#' y <- cbind(FF5F$NA.SMALL.LoBM, FF5F$JP.SMALL.LoBM, FF5F$AP.SMALL.LoBM,
#' FF5F$EU.SMALL.LoBM)
#' ##I fit the data with one bandwidth for each equation
#' ff5f.fit <- tvGLS(x = x, y = y, bw = c(0.89, 1.55, 0.78, 0.31))
#'
#' @method tvGLS list
#' @export

tvGLS.list <- function(x, y, z = NULL, bw, Sigma = NULL, R = NULL, r = NULL,
                       est = c("lc", "ll"), tkernel = c("Epa", "Gaussian"), ...)
{
  if(!is.list(x))
    stop("\n'x' should be a list of matrices. \n")
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  if(tkernel != "Epa" & tkernel != "Gaussian")
    tkernel <- "Epa"
  if(est != "lc" & est != "ll")
    est <- "lc"
  y <- as.matrix(y)
  neq <- ncol(y)
  obs <- nrow(y)
  if (neq != length(x))
    stop("\nIncompatible dimensions. \n")
  for (i in 1:neq)
    x[[i]] <- as.matrix(x[[i]])
  if(is.null(Sigma))
    Sigma <- array(rep(diag(1, neq), obs), dim = c(neq, neq, obs))
  if(length(bw) != 1 & length(bw) != neq)
    stop("\nThere must be a single bandwith or a bandwidth for each equation.\n")
  nvar <- numeric(neq)
  for (i in 1:neq)
    nvar[i] <- ncol(x[[i]])
  y.hat <- matrix(0, obs, neq)
  if(!is.null(z))
    grid <- z
  else
    grid <- (1:obs)/obs
  #nu0 <- integrate(nu, -Inf, Inf,pol = 0,tkernel = tkernel)$value
  if (length(bw)==1)
    bw <- rep(bw, neq)
  resid <- matrix(0, obs, neq)
  theta <- matrix(0, obs, sum(nvar))
  index.length <- numeric(neq)
  for (t in 1:obs)
  {
    x2 <- grid - grid[t]
    y.kernel <- NULL
    x.kernel <- vector ("list", neq)
    eSigma <- eigen(Sigma[,,t], TRUE)
    lambda <- eSigma$values
    if (any(lambda <= 0))
      stop("\n'Sigma' is not positive definite\n")
    A <- diag(lambda^-0.5)  %*%  t(eSigma$vector) %x% Diagonal(obs)
    for (i in 1:neq)
    {
      mykernel <- sqrt(.kernel(x = x2, bw = bw[i], tkernel = tkernel, N = 1))
      y.kernel <- cbind(y.kernel, y[, i]*mykernel)
      xtemp <- x[[i]] * mykernel
      if(est == "ll")
        xtemp <- cbind (xtemp, xtemp * x2)
      x.kernel[[i]] <- xtemp
    }
    x.star <- A %*% Matrix::bdiag(x.kernel)
    y.star <- A %*% as.vector(y.kernel)
    s0 <- Matrix::crossprod(x.star)
    T0 <- Matrix::crossprod(x.star, y.star)
    result <- try(drop(Matrix::solve(s0, T0)), silent = TRUE)
    if(class(result) == "try-error")
      result <- try(drop(Matrix::solve(s0, T0, tol = .Machine$double.xmin)), silent = TRUE)
    if(class(result) == "try-error")
      stop("\nSystem is computationally singular, the inverse cannot be calculated.\n")
    if(!is.null(R))
    {
      Sinv <- Matrix::solve(s0)
      result[1:sum(nvar)] <- result[1:sum(nvar)] - Sinv %*% t(R) %*%
        Matrix::solve(R %*% Sinv %*% t(R)) %*% (R %*% result[1:sum(nvar)]-r)
    }
    theta[t, ] <- result[1:sum(nvar)]
    xt.diag <- as.matrix(Matrix::bdiag(lapply(x,"[",t, ,drop = FALSE)))
    y.hat[t, ] <- xt.diag %*% theta[t, ]
  }
  resid <- as.matrix(y - y.hat)
  return(list( tvcoef = theta, fitted = y.hat, residuals = resid ))
}

#' Estimate Time-varying Coefficients
#'
#' \code{tvGLS} is used to estimate time-varying coefficients SURE using the kernel smoothing GLS
#' @return A list with the estimates, the fitted values and the residuals.
#' @rdname tvGLS
#' @method tvGLS tvsure
#' @export
tvGLS.tvsure <- function(x, ...)
{
  if(class(x) != "tvsure")
    stop ("Function for object of class 'tvsure' \n")
  y <- x$y
  y <- as.matrix(y)
  neq <- ncol(y)
  tkernel <- x$tkernel
  est <- x$est
  if (neq != length(x$x))
    stop("\nIncompatible dimensions. \n")
  if(tkernel != "Epa" & tkernel != "Gaussian")
    tkernel <- "Epa"
  if(est != "lc" & est != "ll")
    est <- "lc"
  obs <- nrow(y)
  bw <- x$bw
  Sigma <- x$Sigma
  nvar <- x$nvar
  est <- x$est
  tkernel <- x$tkernel
  R <- x$R
  r <- x$r
  x <- x$x
  z <- x$z
  for (i in 1:neq)
    x[[i]] <- as.matrix(x[[i]])
  y.hat <- matrix(0, obs, neq)
  if(!is.null(z))
    grid <- z
  else
    grid <- (1:obs)/obs
  if (length(bw) == 1)
    bw <- rep(bw, neq)
  resid <- matrix(0, obs, neq)
  theta <- matrix(0,nrow = obs, ncol = sum(nvar))
  index.length <- numeric(neq)
  bw.N <- bw * obs/2
  kernel.length <- 2 * floor(bw.N)+1
  hmax <- max(kernel.length)
  sqkernel <- matrix(0, nrow = hmax, ncol = neq, byrow = TRUE)
  x2 <- ((hmax:(2 * hmax-1))-(hmax+(2 * hmax-1))*0.5)/obs
  kernel <- list()
  for (i in 1:neq)
  {
    kernel.temp <- .kernel(x2, bw[i], tkernel = tkernel)
    index <- which(kernel.temp!=0)
    sqkernel[index,i] <- sqrt(kernel.temp[index])
  }
  for (t in 1:obs)
  {
    x2 <- grid - grid[t]
    y.kernel <- NULL
    x.kernel <- vector ("list", neq)
    index.kernel <- max(1, floor(max(bw.N)) + 2 - t):min(hmax, obs - t + floor(max(bw.N)) + 1)
    myindex <- max(1,ceiling(t-max(bw.N))):min(obs,floor(t+max(bw.N)))
    eSigma <- eigen(Sigma[,,t], TRUE)
    lambda <- eSigma$values
    if (any(lambda <= 0))
      stop("\n'Sigma' is not positive definite\n")
    A <- diag(lambda^-0.5) %*% t(eSigma$vector) %x% Diagonal(length(index.kernel))
    for (i in 1:neq)
    {
      mykernel <- sqkernel[index.kernel, i]
      y.kernel <- cbind(y.kernel, y[myindex, i]*mykernel)
      xtemp <- x[[i]][myindex, ]*mykernel
      if(est == "ll")
        xtemp <- cbind (xtemp, xtemp * x2[myindex])
      x.kernel[[i]] <- xtemp
    }
    x.star <- A %*% Matrix::bdiag(x.kernel)
    y.star <- A %*% as.vector(y.kernel)
    s0 <- Matrix::crossprod(x.star)
    T0 <- Matrix::crossprod(x.star, y.star)
    result <- try(drop(Matrix::solve(s0, T0)), silent = TRUE)
    if(class(result) == "try-error")
      result <- try(drop(Matrix::solve(s0, T0, tol = .Machine$double.xmin)), silent = TRUE)
    if (!is.null(R))
    {
      Sinv <- Matrix::solve(s0)
      result[1:sum(nvar)] <- result[1:sum(nvar)] - Sinv %*% t(R) %*%
        Matrix::solve(R %*% Sinv %*% t(R)) %*% (R %*% result[1:sum(nvar)]-r)
    }
    theta[t, ] <- result[1:sum(nvar)]
    xt.diag <- as.matrix(Matrix::bdiag(lapply(x,"[",t, ,drop = FALSE)))
    y.hat[t, ] <-  xt.diag %*% theta[t, ]
  }
  resid = y - y.hat
  return(list( tvcoef = theta, fitted = y.hat, residuals = resid ))
}

#' @rdname tvGLS
#' @inheritParams tvGLS
#' @return A list with the estimates, fitted and residuals values.
#' @method tvGLS matrix
#' @export

tvGLS.matrix <- function(x, y, z = NULL, bw, Sigma = NULL, R = NULL, r = NULL, 
                         est = c("lc", "ll"), tkernel = c("Epa", "Gaussian"), ...)
{
  if(!is.matrix(x))
    stop("\n'x' should be a matrix. \n")
  y <- as.numeric(y)
  obs <- nrow(x)
  neq <- length(y)/obs
  if(length(z) != obs)
    stop("\nIncompatible data dimensions.\n")
  if(is.null(Sigma))
    Sigma <- array(rep(diag(1, neq), obs), dim = c(neq, neq, obs))
  if(length(bw) != 1 & length(bw) != neq)
    stop("\nThere must be a single bandwith or a bandwidth for each equation.\n")
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  if(tkernel != "Epa" & tkernel != "Gaussian")
    tkernel <- "Epa"
  if(est != "lc" & est != "ll")
    est <- "lc"
  nvar <- ncol(x)
  y.hat <- matrix(0, obs, neq)
  if(!is.null(z))
    grid <- z
  else
    grid <- (1:obs)/obs
  if (length(bw) == 1)
    bw <- rep(bw, neq)
  resid <- matrix(0, obs, neq)
  theta <- matrix(0, obs, sum(nvar))
  index.length <- numeric(neq)
  bw.N <- bw * obs/2
  kernel.length <- 2 * floor(bw.N)+1
  hmax <- max(kernel.length)
  sqkernel <- matrix(0, nrow = hmax, ncol = neq, byrow = TRUE)
  x2 <- ((hmax:(2 * hmax - 1))-(hmax+(2 * hmax - 1)) * 0.5)/obs

  for (i in 1:neq)
  {
    kernel.temp <- .kernel(x2, bw[i], tkernel = tkernel)
    index <- which(kernel.temp != 0)
    sqkernel[index,i] <- sqrt(kernel.temp[index])
  }
  for (t in 1:obs)
  {
    x2 <- grid -grid[t]
    y.kernel <- NULL
    x.kernel <- vector ("list", neq)
    index.kernel <- max(1, floor(max(bw.N))+2 - t):min(hmax, obs - t + floor(max(bw.N)) + 1)
    myindex <- max(1, ceiling(t - max(bw.N))):min(obs, floor(t + max(bw.N)))
    eSigma <- eigen(Sigma[,,t], TRUE)
    lambda <- eSigma$values
    if (any(lambda <= 0))
      stop("\n'Sigma' is not positive definite\n")
    A <- diag(lambda^-0.5)  %*%  t(eSigma$vector) %x% Diagonal(length(index.kernel))
    for (i in 1:neq)
    {
      mykernel <- sqkernel[index.kernel, i]
      y.kernel <- cbind(y.kernel, y[myindex, i]*mykernel)
      xtemp <- x[[i]][myindex, ]*mykernel
      if(est == "ll")
        xtemp <- cbind (xtemp, xtemp * x2[myindex])
      x.kernel[[i]] <- xtemp
    }
    x.star <- A %*% Matrix::bdiag(x.kernel)
    y.star <- A %*% as.vector(y.kernel)
    s0 <- Matrix::crossprod(x.star)
    T0 <- Matrix::crossprod(x.star, y.star)
    result <- try(drop(Matrix::solve(s0, T0)), silent = TRUE)
    if(class(result) == "try-error")
      result <- try(drop(Matrix::solve(s0, T0, tol = .Machine$double.xmin)), silent = TRUE)
    if (!is.null(R))
    {
      Sinv <- Matrix::solve(s0)
      result <- result[1:sum(nvar)] - Sinv %*% t(R) %*%
        Matrix::solve(R %*% Sinv %*% t(R)) %*% (R %*% result[1:sum(nvar)] - r)
    }
    theta[t, ] <- result[1:sum(nvar)]
    xt.diag <- as.matrix(Matrix::bdiag(lapply(x,"[",t, ,drop = FALSE)))
    y.hat[t, ] <-  xt.diag %*% theta[t, ]
  }
  resid <- as.matrix(y - y.hat)
  return(list( tvcoef = theta, fitted = y.hat, residuals = resid ))
}
