#' Time-Varying Generalised Least Squares
#'
#' \code{tvGLS} estimates time-varying coefficients of SURE using the kernel smoothing GLS.
#'
#'
#' @param x an object used to select a method.
#' @param ... Other arguments passed to specific methods.
#' @return \code{tvGLS} returns a list containing:
#' \item{tvcoef}{An array of dimension obs x nvar x neq (obs = number of observations, nvar = number of variables
#' in each equation, neq = number of equations in the system) with the time-varying coefficients estimates.}
#' \item{fitted}{A matrix of dimension obs x neq with the fited values from the estimation.}
#' \item{residuals}{A matrix of dimension obs x neq with the residuals from the estimation.}
#' @export
#' @import Matrix
#' @import methods

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
#' @param z A vector with the smoothing variable.
#' @param ez (optional) A scalar or vector with the smoothing values. If 
#' values are included then the vector z is used.
#' @param bw A numeric vector.
#' @param Sigma An array.
#' @param R A matrix.
#' @param r A numeric vector.
#' @param est Either "lc" or "ll".
#' @param tkernel Either "Gaussian" or "Epa".
#'
#' @examples
#' data(FF5F)
#' x <- list()
#' ## SMALL/LoBM porfolios time-varying three factor model
#' x[[1]] <- cbind(rep (1, 314), FF5F[, c("NA.Mkt.RF", "NA.SMB",  "NA.HML", "NA.RMW", "NA.CMA")])
#' x[[2]] <- cbind(rep (1, 314), FF5F[, c("JP.Mkt.RF", "JP.SMB",  "JP.HML", "JP.RMW", "JP.CMA")])
#' x[[3]] <- cbind(rep (1, 314), FF5F[, c("AP.Mkt.RF", "AP.SMB",  "AP.HML", "AP.RMW", "AP.CMA")])
#' x[[4]] <- cbind(rep (1, 314), FF5F[, c("EU.Mkt.RF", "EU.SMB",  "EU.HML", "EU.RMW", "EU.CMA")])
#' ##Returns
#' y <- cbind(FF5F$NA.SMALL.LoBM, FF5F$JP.SMALL.LoBM, FF5F$AP.SMALL.LoBM, 
#' FF5F$EU.SMALL.LoBM) 
#' ##Excess returns
#' y <- y - cbind(FF5F$NA.RF, FF5F$JP.RF, FF5F$AP.RF, FF5F$EU.RF)
#' ##I fit the data with one bandwidth for each equation
#' ff5f.fit <- tvGLS(x = x, y = y, bw = c(1.03, 0.44, 0.69, 0.31))
#'
#' @method tvGLS list
#' @export

tvGLS.list <- function(x, y, z = NULL, ez = NULL, bw, Sigma = NULL, R = NULL, r = NULL,
                       est = c("lc", "ll"), tkernel = c("Epa", "Gaussian"), ...)
{
  if(!is.list(x))
    stop("\n'x' should be a list of matrices. \n")
  is.predict <- ifelse (is.null(ez), FALSE, TRUE)
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  if(!(tkernel %in% c("Epa", "Gaussian")))
    tkernel <- "Epa"
  if(!(est %in% c("lc", "ll")))
    est <- "lc"
  y <- as.matrix(y)
  neq <- length(x)
  for (i in 1:neq)
    x[[i]] <- as.matrix(x[[i]])
  obs <- NROW(x[[1]])
  if(!identical(neq, NCOL(y)) | !identical(obs, NROW(y)))
    stop("The number of equations in 'x' and 'y' are different \n")
  if(is.null(Sigma))
    Sigma <- array(rep(diag(1, neq), obs), dim = c(neq, neq, obs))
  if(length(bw) != 1 & length(bw) != neq)
    stop("\nThere must be a single bandwith or a bandwidth for each equation.\n")
  nvar <- numeric(neq)
  for (i in 1:neq)
    nvar[i] <- NCOL(x[[i]])
  if(!is.null(z))
  {
    if(length(z) != obs)
    {
      stop("\nDimensions of 'y' and 'z' are not compatible.\n")
    }
    grid <- z
  }
  else
    grid <- (1:obs)/obs
  if (is.null(ez))
    ez <- grid
  eobs <- NROW(ez)
  if (length(bw) == 1)
    bw <- rep(bw, neq)
  if(!is.null(R))
  {
    R <- as.matrix(R)
    if(NCOL(R) != sum(nvar))
      stop("\nWrong dimension of 'R', it should have as many columns as variables 
           in the whole system. \n")
    if (is.null(r))
      r <- rep(0, NROW(R))
    else if (length(r) == 1)
      r <- rep(r, NROW(R))
    else if (length(r) != NROW(R) & length(r) != 1)
      stop("\nWrong dimension of 'r', it should be as long as the number of 
           rows in 'R'. \n")
  }
  y.hat <- matrix(0, eobs, neq)
  resid <- matrix(0, eobs, neq)
  theta <- matrix(0, eobs, sum(nvar))
  for (t in 1:eobs)
  {
    x2 <- grid - ez[t]
    y.kernel <- NULL
    x.kernel <- vector ("list", neq)
    eSigma <- eigen(Sigma[,,t], TRUE)
    lambda <- eSigma$values
    if (any(lambda <= 0))
      stop("\n'Sigma' is not positive definite.\n")
    A <- diag(lambda^-0.5) %*% t(eSigma$vectors) %x% Matrix::Diagonal(obs)
    for (i in 1:neq)
    {
      mykernel <- sqrt(.kernel(x = x2, bw = bw[i], tkernel = tkernel, N = 1))
      y.kernel <- cbind(y.kernel, y[, i] * mykernel)
      xtemp <- x[[i]] * mykernel
      if(est == "ll")
        xtemp <- cbind(xtemp, xtemp * x2)
      x.kernel[[i]] <- xtemp
    }
    x.star <- A %*% Matrix::bdiag(x.kernel)
    y.star <- A %*% methods::as(y.kernel, "sparseVector")
    s0 <- Matrix::crossprod(x.star)
    T0 <- Matrix::crossprod(x.star, y.star)
    result <- try(drop(Matrix::solve(s0, T0)), silent = TRUE)
    if(class(result) == "try-error")
      result <- try(drop(Matrix::solve(s0, T0, tol = .Machine$double.xmin)), silent = TRUE)
    if(class(result) == "try-error")
      stop("\nSystem is computationally singular, the inverse cannot be calculated. 
           Possibly, the 'bw' is too small for values in 'ez'.\n")
    if(est =="ll")
    {
      temp <- result[1:nvar[1]]
      for(i in 2:neq)
        temp <- c(temp, result[(1:nvar[i]) + 2 * sum(nvar[1:(i-1)])])
      theta[t, ] <- temp
    }
    else
      theta[t, ] <- result
    if(!is.null(R))
    {
      Sinv <- Matrix::solve(s0)
      theta[t, ] <- theta[t, ] - Sinv %*% t(R) %*%
        Matrix::solve(R %*% Sinv %*% t(R)) %*% (R %*% theta[t, ] - r)
    }
    xt.diag <- as.matrix(Matrix::bdiag(lapply(x,"[",t, ,drop = FALSE)))
    y.hat[t, ] <- xt.diag %*% theta[t, ]
  }
  if(!is.predict)
    resid <- as.matrix(y - y.hat)
  return(list( tvcoef = theta, fitted = y.hat, residuals = resid ))
}

#' @rdname tvGLS
#' @inheritParams tvGLS
#' @method tvGLS tvsure
#' @export
tvGLS.tvsure <- function(x, ...)
{
  if(!inherits(x, c("tvsure")))
    stop ("Function for object of class 'tvsure' \n")
  y <- x$y
  z <- x$z
  ez <- x$ez
  bw <- x$bw
  Sigma <- x$Sigma
  neq <- x$neq
  R <- x$R
  r <- x$r
  est <- x$est
  tkernel <- x$tkernel
  x <- x$x
  result <- tvGLS(x, y, z, ez, bw, Sigma, R, r, est, tkernel)
  return(result)
}

#' @rdname tvGLS
#' @inheritParams tvGLS
#' @method tvGLS matrix
#' @export

tvGLS.matrix <- function(x, y, z = NULL, ez = NULL, bw, Sigma = NULL, 
                         R = NULL, r = NULL, est = c("lc", "ll"), 
                         tkernel = c("Epa", "Gaussian"), ...)
{
  if(!is.matrix(x))
    stop("\n'x' should be a matrix. \n")
  is.predict <- ifelse (is.null(ez), FALSE, TRUE)
  y <- as.numeric(y)
  obs <- NROW(x)
  neq <- length(y)/obs
  if (neq != NCOL(x))
    stop("\nNumber of equations in 'x' and 'y' are not compatible.\n")
  if(is.null(Sigma))
    Sigma <- array(rep(diag(1, neq), obs), dim = c(neq, neq, obs))
  if(length(bw) != 1 & length(bw) != neq)
    stop("\nThere must be a single bandwith or a bandwidth for each equation.\n")
  if (length(bw) == 1)
    bw <- rep(bw, neq)
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  if(!(tkernel %in% c("Epa", "Gaussian")))
    tkernel <- "Epa"
  if(!(est %in% c("lc", "ll")))
    est <- "lc"
  nvar <- NCOL(x)
  if(!is.null(z))
  {
    if(length(z) != obs)
      stop("\nDimensions of 'x' and 'z' are not compatible.\n")
    grid <- z
  }
  else
    grid <- (1:obs)/obs
  if (is.null(ez))
    ez <- grid
  eobs <- NROW(ez)
  y.hat <- matrix(0, eobs, neq)
  resid <- matrix(0, eobs, neq)
  theta <- matrix(0, eobs, sum(nvar))
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
  for (t in 1:eobs)
  {
    x2 <- grid - ez[t]
    y.kernel <- NULL
    x.kernel <- vector ("list", neq)
    index.kernel <- max(1, floor(max(bw.N))+2 - t):min(hmax, obs - t + floor(max(bw.N)) + 1)
    myindex <- max(1, ceiling(t - max(bw.N))):min(obs, floor(t + max(bw.N)))
    eSigma <- eigen(Sigma[,,t], TRUE)
    lambda <- eSigma$values
    if (any(lambda <= 0))
      stop("\n'Sigma' is not positive definite\n")
    A <- diag(lambda^-0.5)  %*%  t(eSigma$vector) %x% Matrix::Diagonal(length(index.kernel))
    for (i in 1:neq)
    {
      mykernel <- sqkernel[index.kernel, i]
      y.kernel <- cbind(y.kernel, y[myindex, i]*mykernel)
      xtemp <- x[[i]][myindex, ]*mykernel
      if(est == "ll")
        xtemp <- cbind (xtemp, xtemp * x2[myindex])
      x.kernel[[i]] <- xtemp
    }
    x.star <- try(A %*% Matrix::bdiag(x.kernel), silent = TRUE)
    y.star <- try(A %*% methods::as(y.kernel, "sparseVector"), silent = TRUE)
    s0 <- Matrix::crossprod(x.star)
    T0 <- Matrix::crossprod(x.star, y.star)
    result <- try(drop(Matrix::solve(s0, T0)), silent = TRUE)
    if(class(result) == "try-error")
      result <- try(drop(Matrix::solve(s0, T0, tol = .Machine$double.xmin)), silent = TRUE)
    if(class(result) == "try-error")
      stop("\nSystem is computationally singular, the inverse cannot be calculated. 
           Possibly, the 'bw' is too small for values in 'ez'.\n")
    if(est =="ll")
    {
      temp <- result[1:nvar[1]]
      for(i in 2:neq)
        temp <- c(temp, result[(1:nvar[i]) + 2 * sum(nvar[1:(i-1)])])
      theta[t, ] <- temp
    }
    else
      theta[t, ] <- result
    if(!is.null(R))
    {
      Sinv <- Matrix::solve(s0)
      theta[t, ] <- theta[t, ] - Sinv %*% t(R) %*%
        Matrix::solve(R %*% Sinv %*% t(R)) %*% (R %*% theta[t, ] - r)
    }
    xt.diag <- as.matrix(Matrix::bdiag(lapply(x,"[",t, ,drop = FALSE)))
    y.hat[t, ] <-  xt.diag %*% theta[t, ]
  }
  if(!is.predict)
    resid <- as.matrix(y - y.hat)
  return(list( tvcoef = theta, fitted = y.hat, residuals = resid ))
}
