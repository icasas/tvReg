#' Time-Varying Random Effects Estimation
#'
#' \code{tvRE} estimate time-varying coefficient of a random effects 
#' panel data model using kernel smoothing.
#'
#' @param x An object used to select a method.
#' @param ... Other arguments passed to specific methods.
#' @return \code{tvRE} returns a list containing:
#' \item{coefficients}{A vector of length obs, number of observations with
#' the time-varying estimates.}
#' \item{fitted}{A vector of length obs with the fited values from the estimation.}
#' \item{residuals}{A vector of length obs with the residuals from the estimation.}
#' \item{alpha}{A vector of length neq with the fixed effects.}
#' @export
tvRE <- function(x, ...) UseMethod("tvRE", x)

#' @rdname tvRE
#' @method tvRE matrix
#' @param y A vector with dependent variable.
#' @param z A vector with the variable over which coefficients are smooth over.
#' @param ez (optional) A scalar or vector with the smoothing values. If 
#' values are included then the vector z is used.
#' @param bw A numeric vector with the bandwidth.
#' @param Sigma NULL (default) or a matrix of size obs x obs..
#' @param neq A scalar with the number of equations
#' @param obs A scalar with the number of time observations
#' @param est The nonparametric estimation method, one of "lc" (default) for linear constant
#'  or "ll" for local linear.
#' @param tkernel A character, either "Triweight" (default), "Epa" or "Gaussian" kernel function.
#
#' @export
tvRE.matrix <-function(x, y, z = NULL, ez = NULL, bw, Sigma = NULL, neq, obs,
               est = c("lc", "ll"), tkernel = c("Triweight", "Epa", "Gaussian"), ...)
{
  x <- as.matrix(x)
  y <- as.numeric(y)
  if(!identical(NROW(y), NROW(x)))
    stop("\nDimensions of 'x' and 'y' are not compatible.\n")
  if(!is.numeric(bw))
    stop ("Argument 'bw' should be a scalar. \n")
  if(is.null(Sigma))
    Sigma <- diag(1, obs)
  is.predict <- ifelse (is.null(ez), FALSE, TRUE)
  if(!is.null(z))
  {
    if(length(z) != obs)
      stop("\nDimensions of 'x' and 'z' are not compatible. \nThe size of 'z' should be 'obs'.\n")
    grid <- z
  }
  else
    grid <- (1:obs)/obs
  if (is.null(ez))
    ez <- grid
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  nvar <- NCOL(x)
  eobs <- NROW(ez)
  fitted = resid = alpha <- numeric(neq*eobs)
  theta <- matrix(0, eobs, nvar)
  alpha.i <- numeric(neq)
  xnew <- matrix(NA, neq*obs, nvar)
  if(est == "ll")
    xnew <- matrix(NA, neq*obs, 2*nvar)
  ynew <- y
  for (t in 1:eobs)
  {
    tau0 <- grid - ez[t]
    kernel.bw <- sqrt(.kernel(tau0, bw, tkernel))
    if (length (kernel.bw != 0) < 3)
      stop("Bandwidth is too small.\n")
    xtemp <- x
    if (est == "ll")
      xtemp <- cbind(xtemp, xtemp * tau0)
    eSigma <- eigen(Sigma, TRUE)
    if (any(eSigma$value <= 0))
      stop("\n'Sigma' is not positive definite.\n")
    A <- diag(eSigma$values^-0.5) %*% t(eSigma$vectors)
    for (i in 1:neq)
    {
      ind <- (i-1)*obs + (1:obs)
      xnew[ind, ] <- A%*%(xtemp[ind,] * kernel.bw)
      ynew[ind] <- A%*%(y[ind] * kernel.bw)
    }
    result <- try(qr.solve(crossprod(xnew), crossprod(xnew, ynew)), silent =TRUE)
    if(inherits(result, "try-error"))
      result <- try(qr.solve(crossprod(xnew), crossprod(xnew, ynew), tol = .Machine$double.xmin), 
                    silent = TRUE)
    if(inherits(result, "try-error"))
      stop("\nSystem is computationally singular, the inverse cannot be calculated. 
           Possibly, the 'bw' is too small for values in 'ez'.\n")
    theta[t,] <- result [1:nvar]
    for (i in 1:neq)
      fitted[t+(i-1)*eobs] <- sum(x[t+(i-1)*obs, ]*theta[t,])
    alpha[t+((1:neq)-1)*eobs] <- y[t+((1:neq)-1)*obs] - fitted[t+((1:neq)-1)*eobs]
  }
  for(i in 1:neq)
    alpha.i[i] <- mean(alpha[(i-1)*obs + 1:obs])
  if(!is.predict)
  {
    resid <- y - fitted - rep(alpha.i, each = obs)
    resid.mat <- matrix(resid^2, nrow = obs, ncol = neq, byrow = TRUE)
    u.var <- sum(apply((resid.mat - apply(resid.mat, 2, mean))^2, 2, sum))/(neq*obs - neq)
    Sigma <- matrix(u.var, obs, obs)
    diag(Sigma) <- mean(apply(resid.mat, 1, mean))
  }
  else
    fitted = resid <- NULL
  return(list(coefficients = theta, fitted = fitted, residuals = resid, 
              alpha = alpha.i, Sigma = Sigma))
}

#' @rdname tvRE
#' @method tvRE tvplm
#' @export
tvRE.tvplm <- function(x, ...)
{
  return(tvRE(x = x$x, x$y, x$z, x$ez, x$bw, x$Sigma, x$neq, x$obs,
              x$est, x$tkernel))
}


#' @rdname tvReg-internals
#' @keywords internal
.tvRE.cv<-function(bw, x, y, z = NULL, neq, obs, cv.block = 0,
                   est = c("lc", "ll"), 
                   tkernel = c("Triweight", "Epa", "Gaussian"))
{
  x <- as.matrix(x)
  nvar <- NCOL (x)
  fitted <- numeric(neq*obs)
  if(!is.null(z))
  {
    if(length(z) != obs)
      stop("\nDimensions of 'x' and 'z' are not compatible. \nThe size of 'z' should be 'obs'.\n")
    grid <- z
  }
  else
    grid <- (1:obs)/obs
  xnew <- matrix(NA, neq*obs, nvar)
  if(est == "ll")
    xnew <- matrix(NA, neq*obs, 2*nvar)
  ynew <- y
  for (t in 1:obs)
  {
    tau0 <- grid - grid[t]
    kernel.bw <- sqrt(.kernel(x = tau0, bw = bw, tkernel = tkernel))
    kernel.bw[max(1, (t- cv.block)):min((t+cv.block), obs)] <- 0
    if (sum(kernel.bw != 0) < 3)
      return (.Machine$double.xmax)
    xtemp <- x
    if (est=="ll")
      xtemp <- cbind(xtemp, xtemp * tau0)
    for (i in 1:neq)
    {
      ind <- (i-1)*obs + (1:obs)
      xnew[ind, ] <- xtemp[ind,] * kernel.bw
      ynew[ind] <- y[ind] * kernel.bw
    }
    theta <- (solve(crossprod(xnew))%*%crossprod(xnew, ynew))[1:nvar,]
    for (i in 1:neq)
      fitted[t+(i-1)*obs] <- sum(x[t+(i-1)*obs,]*theta)
  }
  resid <- y - fitted
  return(mean(resid^2))
}

