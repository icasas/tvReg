#' Time-Varying Ordinary Least Squares
#'
#' \code{tvOLS} estimate time-varying coefficient of univariate 
#' linear models using the kernel smoothing OLS.
#'
#' @param x an object used to select a method.
#' @param ... Other arguments passed to specific methods.
#' @return \code{tvGLS} returns a list containing:
#' \item{tvcoef}{A vector of length obs, number of observations with
#' the time-varying estimates.}
#' \item{fitted}{A vector of length obs with the fited values from the estimation.}
#' \item{residuals}{A vector of length obs with the residuals from the estimation.}
#' @export
#' @import methods
tvOLS <- function(x, ...) UseMethod("tvOLS", x)

#' @param x A matrix with all regressors.
#' @param y A vector with dependent variable.
#' @param z A vector with the variable over which coefficients are smooth over.
#' @param ez (optional) A scalar or vector with the smoothing values. If 
#' values are included then the vector z is used.
#' @param bw A numeric vector.
#' @inheritParams tvSURE
#' @param singular.ok	Logical. If FALSE, a singular model is an error.
#' @examples
#' tau <- seq(1:500)/500
#' beta <- data.frame(beta1 = sin(2*pi*tau), beta2= 2*tau)
#' X <- data.frame(X1 = rnorm(500), X2 = rchisq(500, df = 4))
#' error <- rt(500, df = 10)
#' y <- apply(X*beta, 1, sum) + error
#' coef.lm <- stats::lm(y~0+X1+X2, data = X)$coef
#' coef.tvlm <-  tvOLS(x = X, y = y, bw = 0.1)$tvcoef
#' plot(tau,beta[, 1], type="l", main="", ylab = expression(beta[1]), xlab = expression(tau),
#' ylim = range(beta[,1], coef.tvlm[, 1]))
#' abline(h = coef.lm[1], col = 2)
#' lines(tau, coef.tvlm[, 1], col = 4)
#' legend("topright", c(expression(beta[1]), "lm", "tvlm"), col = c(1, 2, 4), bty="n", lty = 1)
#'
#' @seealso \code{\link{bw}} for bandwidth selection, \code{\link{tvLM}} and
#' \code{\link{tvAR}}.
#' @method tvOLS matrix
#' @export

tvOLS.matrix <- function(x, y, z = NULL, ez = NULL, bw, est = c("lc", "ll"), 
                  tkernel = c("Epa", "Gaussian"), 
                  singular.ok = singular.ok, ...)
{
  x <- as.matrix(x)
  y <- as.numeric(y)
  is.predict <- ifelse (is.null(ez), FALSE, TRUE)
  obs <- NROW(x)
  if(!identical(length(y), obs))
    stop("\nDimensions of 'x' and 'y' are not compatible.\n")
  if(!is.numeric(bw))
    stop ("Argument 'bw' should be a scalar. \n")
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
      stop("\nDimensions of 'x' and 'z' are not compatible\n")
    grid <- z
  }
  else
    grid <- (1:obs)/obs
  if (is.null(ez))
    ez <- grid
  eobs <- NROW(ez)
  fitted <- numeric(eobs)
  resid <- numeric(eobs)
  theta <- matrix(0, eobs, nvar)
  for (t in 1:eobs)
  { 
    x2 <- grid - ez[t]
    kernel.bw <- .kernel(x = x2, bw = bw, tkernel = tkernel, N = 1)
    myindex <- which(kernel.bw != 0)
    if (length (myindex) < 3)
      stop("Bandwidth is too small or 'ez' is too large.\n")
    if (sum(myindex) == 0)
      return (.Machine$double.xmax)
    xtemp <- x[myindex, ]
    if (est=="ll")
      xtemp <- cbind(xtemp, xtemp * x2[myindex])
    ytemp <- y[myindex]
    result <- stats::lm.wfit(x = as.matrix(xtemp), y = ytemp, w = kernel.bw[myindex])
    theta[t,] <- result$coef[1:nvar]
    fitted[t] <- crossprod(x[t, !is.na(theta[t,])], theta[t, !is.na(theta[t,])])
  }
  if(!is.predict)
    resid <- y - fitted
  return(list(tvcoef = theta, fitted = fitted, residuals = resid))
}

#' @rdname tvOLS
#' @method tvOLS tvlm
#' @export
tvOLS.tvlm <- function(x, ...)
{
  if(!inherits(x, c("tvlm", "tvar")))
    stop ("Function for object of class 'tvlm' and 'tvar'. \n")
  y <- x$y
  y <- as.matrix(y)
  tkernel <- x$tkernel
  est <- x$est
  obs <- x$obs
  bw <- x$bw
  z <- x$z
  ez <- x$ez
  singular.ok <- x$singular.ok
  return(tvOLS(x = x$x, y, z, ez, bw, est, tkernel, singular.ok))
}
  
#' @rdname tvOLS
#' @method tvOLS tvar
#' @export
tvOLS.tvar <- tvOLS.tvlm

#' @rdname tvOLS
#' @method tvOLS tvvar
#' @export
tvOLS.tvvar <- function(x, ...)
{
  if(!inherits(x, c("tvvar")))
    stop ("Function for object of class 'tvvar'. \n")
  neq <- x$neq
  yend <- x$y
  rhs <- x$x
  z <- x$z
  ez <- x$ez
  bw <- x$bw
  est <- x$est
  tkernel <- x$tkernel
  singular.ok <- x$singular.ok
  equation <- list()
  resid = fitted <- matrix(0, nrow = x$obs, ncol = neq)
  for (i in 1:neq)
  {
    y <- yend[, i]
    results <- tvOLS(x = rhs, y = y, z = z, ez = ez, bw = bw[i], 
                     est = est, tkernel = tkernel, singular.ok = singular.ok)
    equation[[colnames(yend)[i]]] <- results$tvcoef
    colnames(equation[[colnames(yend)[i]]]) <- colnames(rhs)
    resid[,i] <- results$residuals
    fitted[, i] <- results$fitted
  }
  
  return(list(tvcoef = equation, fitted = fitted, residuals = resid))
}


#' @rdname tvReg-internals
#' @keywords internal
.tvOLS.cv <- function(bw, x, y, z = NULL, est = c("lc", "ll"), 
                      tkernel = c("Epa", "Gaussian"), singular.ok = TRUE)
{
  x <- as.matrix(x)
  obs <- NROW(x)
  if(!identical(length(y), obs))
    stop("\nDimensions of 'x' and 'y' are not compatible\n")
  fitted <- resid.2 <- numeric(obs)
  nvar <- NCOL(x)
  bw <- abs(bw)
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  if(!(tkernel %in% c("Epa", "Gaussian")))
    tkernel <- "Epa"
  if(!(est %in% c("lc", "ll")))
    est <- "lc"
  if(!is.null(z))
  {
    if(length(z) != obs)
      stop("\nDimensions of 'x' and 'z' are not compatible\n")
    grid <- z
  }
  else
    grid <- (1:obs)/obs
  theta <- matrix(0, obs, nvar)
  fitted <- numeric(obs)
  for (t in 1:obs)
  {
    x2 <- grid - grid[t]
    kernel.bw <- .kernel(x = x2, bw = bw, tkernel = tkernel, N = 1)
    kernel.bw[t] <- 0
    myindex <- which(kernel.bw != 0)
    if (sum(myindex != 0) < 3)
      return (.Machine$double.xmax)
    xtemp <- x[myindex, ]
    if (est == "ll")
      xtemp <- cbind(xtemp, xtemp * x2[myindex])
    ytemp <- y[myindex]
    result <- stats::lm.wfit(x = as.matrix(xtemp), y = ytemp, w = kernel.bw[myindex],
                             singular.ok = singular.ok)
    theta[t,] <- result$coef[1:nvar]
    fitted[t] <- crossprod(x[t, !is.na(theta[t,])], theta[t, !is.na(theta[t,])])
  }
  mse <- mean((y - fitted)^2)
  return(mse)
}
