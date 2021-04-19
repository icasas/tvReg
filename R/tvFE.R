
#' Time-Varying Fixed Effects Estimation
#'
#' \code{tvFE} estimate time-varying coefficient of fixed effects 
#' panel data models using kernel smoothing.
#'
#' @param x An object used to select a method.
#' @param ... Other arguments passed to specific methods.
#' @return \code{tvFE} returns a list containing:
#' \item{coefficients}{A vector of length obs, number of observations with
#' the time-varying estimates.}
#' \item{fitted}{A vector of length obs with the fited values from the estimation.}
#' \item{residuals}{A vector of length obs with the residuals from the estimation.}
#' \item{alpha}{A vector of length neq with the fixed effects.}
#' @export
#' @import methods
tvFE <- function(x, ...) UseMethod("tvFE", x)

#' @rdname tvFE
#' @method tvFE matrix
#' @param y A vector with dependent variable.
#' @param z A vector with the variable over which coefficients are smooth over.
#' @param ez (optional) A scalar or vector with the smoothing values. If 
#' values are included then the vector z is used.
#' @param bw A numeric vector.
#' @param neq A sclar with the number of equations
#' @param obs A scalar with the number of time observations
#' @param est The nonparametric estimation method, one of "lc" (default) for linear constant
#'  or "ll" for local linear.
#' @param tkernel A character, either "Triweight" (default), "Epa" or "Gaussian" kernel function.
#
#' @export
tvFE.matrix<-function(x, y, z = NULL, ez = NULL, bw, neq, obs,
               est = c("lc", "ll"), tkernel = c("Triweight", "Epa", "Gaussian"), ...)
{
  x <- as.matrix(x)
  y <- as.numeric(y)
  if(!identical(NROW(y), NROW(x)))
    stop("\nDimensions of 'x' and 'y' are not compatible.\n")
  if(!is.numeric(bw))
    stop ("Argument 'bw' should be a scalar. \n")
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
  fitted = resid <- numeric(neq*eobs)
  theta <- matrix(0, eobs, nvar)
  alpha <- matrix(0, nrow = eobs, ncol = neq-1)
  for (t in 1:eobs)
  {
    tau0 <- grid - ez[t]
    kernel.bw <- .kernel(tau0, bw, tkernel)
    if (length (kernel.bw != 0) < 3)
      stop("Bandwidth is too small.\n")
    D <- t(cbind(rep(-1, neq-1), diag(1, neq-1)))%x%rep(1, obs)
    xtemp <- x
    if(est == "ll")
      xtemp <- cbind (xtemp, xtemp * tau0)
    WH <- diag(neq)%x%diag(kernel.bw)
    DW <- crossprod(D, WH)
    temp <- qr.solve(DW%*%D)%*%DW
    MH <- diag(neq*obs) - D %*% temp
    SH <- crossprod(MH, WH)%*%MH
    x.tilde <- crossprod(xtemp, SH)
    T0 <- x.tilde %*% y
    s0 <- x.tilde %*% xtemp
    result <- try(qr.solve(s0, T0), silent = TRUE)
    if(inherits(result, "try-error"))
      result <- try(qr.solve(s0, T0, tol = .Machine$double.xmin ), silent = TRUE)
    if(inherits(result, "try-error"))
      stop("\nSystem is computationally singular, the inverse cannot be calculated. 
           Possibly, the 'bw' is too small for values in 'ez'.\n")
    theta[t, ] <- result[1:nvar]
    xtheta <- x%*%theta[t,]
    alpha[t,] <- as.numeric(temp%*%(y - xtheta))
    fitted[t+((1:neq)-1)*eobs] <- xtheta[t+((1:neq)-1)*obs]
  }
  alpha <- apply(alpha, 2, mean)
  if(!is.predict)
  {
    resid <- y - fitted 
    D <- Matrix(t(cbind(rep(-1, neq-1), diag(1, neq-1)))%x%rep(1, eobs))
    fitted <- drop(fitted + D %*% alpha)
  }
  else
    fitted = resid <- NULL
    
  return(list(coefficients = theta, fitted = fitted, residuals = resid, alpha = c(-sum(alpha), alpha)))
}

#' @rdname tvFE
#' @method tvFE tvplm
#' @export
tvFE.tvplm <- function(x, ...)
{
  return(tvFE(x = x$x, x$y, x$z, x$ez, x$bw, x$neq, x$obs,
                x$est, x$tkernel))
}

#' @rdname tvReg-internals
#' @keywords internal
.tvFE.cv<-function(bw, x, y, z = NULL, neq, obs, cv.block = 0,
                   est = c("lc", "ll"), 
                   tkernel = c("Triweight", "Epa", "Gaussian"))
{
  x <- as.matrix(x)
  fitted = resid <- numeric(obs*neq)
  nvar <- NCOL(x)
  alpha <- matrix(0, nrow = obs, ncol = neq-1) 
  if(!is.null(z))
  {
    if(length(z) != obs)
      stop("\nDimensions of 'x' and 'z' are not compatible. \nThe size of 'z' should be 'obs'.\n")
    grid <- z
  }
  else
    grid <- (1:obs)/obs
  for (t in 1:obs)
  {
    tau0 <- grid - grid[t]
    kernel.bw <- .kernel(tau0, bw, tkernel)
    kernel.bw[max(1, t-cv.block):min(t+cv.block, obs)] <- 0
    if (sum(kernel.bw != 0) < 3) 
      return (.Machine$double.xmax)
    D <- t(cbind(rep(-1, neq-1), diag(1, neq-1)))%x%rep(1, obs)
    xtemp <- x
    if(est == "ll")
      xtemp <- cbind (xtemp, xtemp * tau0)
    WH <- diag(neq)%x%diag(kernel.bw)
    DW <- crossprod(D, WH)
    temp <- qr.solve(DW%*%D)%*%DW
    MH <- diag(neq*obs) - D %*% temp
    SH <- crossprod(MH, WH)%*%MH
    x.tilde <- crossprod(xtemp, SH)
    T0 <- x.tilde %*% y
    s0 <- x.tilde %*% xtemp
    result <- try(qr.solve(s0, T0), silent = TRUE)
    if(inherits(result, "try-error"))
      return (.Machine$double.xmax)
    xtheta <- x%*%matrix(result[1:nvar])
    alpha[t,] <- as.numeric(temp%*%(y - xtheta))
    fitted[t+((1:neq)-1)*obs] <- xtheta[t+((1:neq)-1)*obs]
  }
  alpha <- apply(alpha, 2, mean)
  D <- t(cbind(rep(-1, neq-1), diag(1, neq-1)))%x%rep(1, obs)
  resid <- y - fitted - D%*%alpha
  return(mean(resid^2))
}

