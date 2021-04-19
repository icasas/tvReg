#' Time-Varying Ordinary Least Squares
#'
#' \code{tvOLS} estimate time-varying coefficient of univariate 
#' linear models using the kernel smoothing OLS.
#'
#' @param x An object used to select a method.
#' @param ... Other arguments passed to specific methods.
#' @return \code{tvOLS} returns a list containing:
#' \item{coefficients}{A vector of length obs, number of observations time observations.}
#' \item{fitted}{A vector of length obs with the fitted values from the estimation.}
#' \item{residuals}{A vector of length obs with the residuals from the estimation.}
#' @export
#' @import methods
tvOLS <- function(x, ...) UseMethod("tvOLS", x)

#' @rdname tvOLS
#' @method tvOLS matrix
#' @param y A vector with dependent variable.
#' @param z A vector with the variable over which coefficients are smooth over.
#' @param ez (optional) A scalar or vector with the smoothing values. If 
#' values are included then the vector z is used.
#' @param bw A numeric vector.
#' @inheritParams tvLM
#' @param singular.ok	Logical. If FALSE, a singular model is an error.
#' @examples
#' tau <- seq(1:500)/500
#' beta <- data.frame(beta1 = sin(2*pi*tau), beta2 = 2*tau)
#' X <- data.frame(X1 = rnorm(500), X2 = rchisq(500, df = 4))
#' error <- rt(500, df = 10)
#' y <- apply(X*beta, 1, sum) + error
#' coef.lm <- stats::lm(y~0+X1+X2, data = X)$coef
#' coef.tvlm <-  tvOLS(x = as.matrix(X), y = y, bw = 0.1)$coefficients
#' plot(tau, beta[, 1], type="l", main="", ylab = expression(beta[1]), xlab = expression(tau),
#' ylim = range(beta[,1], coef.tvlm[, 1]))
#' abline(h = coef.lm[1], col = 2)
#' lines(tau, coef.tvlm[, 1], col = 4)
#' legend("topright", c(expression(beta[1]), "lm", "tvlm"), col = c(1, 2, 4), bty="n", lty = 1)
#'
#' @seealso \code{\link{bw}} for bandwidth selection, \code{\link{tvLM}} and
#' \code{\link{tvAR}}.
#' @export
tvOLS.matrix <- function(x, y, z = NULL, ez = NULL, bw, est = c("lc", "ll"), 
                  tkernel = c("Triweight", "Epa", "Gaussian"), singular.ok = TRUE, ...)
{
  x <- as.matrix(x)
  y <- as.numeric(y)
  obs <- NROW(x)
  if(!identical(length(y), obs))
    stop("\nDimensions of 'x' and 'y' are not compatible.\n")
  if(!is.numeric(bw))
    stop ("Argument 'bw' should be a scalar. \n")
  is.predict <- ifelse (is.null(ez), FALSE, TRUE)
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
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  nvar <- NCOL(x)
  eobs <- NROW(ez)
  fitted = resid <- numeric(eobs)
  theta <- matrix(0, eobs, nvar)
  for (t in 1:eobs)
  { 
    tau0 <- grid - ez[t]
    kernel.bw <- .kernel(x = tau0, bw = bw, tkernel = tkernel)
    k.index <- which(kernel.bw != 0)
    if (length (k.index) < 3)
      stop("Bandwidth is too small for values in 'ez'.\n")
    xtemp <- x[k.index, ]
    if (est=="ll")
      xtemp <- cbind(xtemp, xtemp * tau0[k.index])
    result <- stats::lm.wfit(x = as.matrix(xtemp), y = y[k.index], 
                             w = kernel.bw[k.index], singular.ok = singular.ok)
    theta[t,] <- result$coefficients[1:nvar]
    fitted[t] <- crossprod(x[t, !is.na(theta[t,])], theta[t, !is.na(theta[t,])])
  }
  if(!is.predict)
    resid <- y - fitted
  else
    fitted = resid <- NULL
  return(list(coefficients = theta, fitted = fitted, residuals = resid))
}

#' @rdname tvOLS
#' @method tvOLS tvlm
#' @export
tvOLS.tvlm <- function(x, ...)
{
  return(tvOLS(x = x$x, x$y, x$z, x$ez, x$bw, x$est, x$tkernel, x$singular.ok))
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
  equation <- list()
  neq <- x$neq
  rhs <- x$x
  eqnames <- colnames(x$y)
  is.predict <- ifelse (is.null(x$ez), FALSE, TRUE)
  resid = fitted <- NULL
  if(!is.predict)
    resid = fitted <- matrix(0, nrow = x$obs, ncol = neq)
  for (i in 1:neq)
  {
    results <- tvOLS(x = rhs, y = x$y[, i], z = x$z, ez = x$ez, bw = x$bw[i], est = x$est, tkernel = x$tkernel, 
                     singular.ok = x$singular.ok)
    equation[[eqnames[i]]] <- results$coefficients
    colnames(equation[[eqnames[i]]]) <- colnames(rhs)
    if(!is.predict)
    {
      resid[, i] <- results$residuals
      fitted[, i] <- results$fitted
    }
  }
  return(list(coefficients = equation, fitted = fitted, residuals = resid))
}


#' @rdname tvReg-internals
#' @keywords internal
.tvOLS.cv <- function(bw, x, y, z = NULL,  cv.block = 0, est = c("lc", "ll"), 
                      tkernel = c("Triweight", "Epa", "Gaussian"), singular.ok = TRUE)
{
  x <- as.matrix(x)
  obs <- NROW(x)
  if(!identical(length(y), obs))
    stop("\nDimensions of 'x' and 'y' are not compatible\n")
  if(!is.null(z))
  {
    if(length(z) != obs)
      stop("\nDimensions of 'x' and 'z' are not compatible\n")
    grid <- z
  }
  else
    grid <- (1:obs)/obs
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  fitted = resid.2 <- numeric(obs)
  nvar <- NCOL(x)
  bw <- abs(bw)
  theta <- matrix(0, obs, nvar)
  fitted <- numeric(obs)
  for (t in 1:obs)
  {
    tau0 <- grid - grid[t]
    kernel.bw <- .kernel(x = tau0, bw = bw, tkernel = tkernel)
    kernel.bw[max(1, t - cv.block):min(t + cv.block, obs)] <- 0
    k.index <- which(kernel.bw != 0)
    if (sum(k.index != 0) < 3)
      return (.Machine$double.xmax)
    xtemp <- x[k.index, ]
    if (est == "ll")
      xtemp <- cbind(xtemp, xtemp * tau0[k.index])
    result <- stats::lm.wfit(x = as.matrix(xtemp), y = y[k.index], w = kernel.bw[k.index],
                             singular.ok = singular.ok)
    theta[t,] <- result$coefficients[1:nvar]
    fitted[t] <- crossprod(x[t, !is.na(theta[t,])], theta[t, !is.na(theta[t,])])
  }
  mse <- mean((y - fitted)^2)
  return(mse)
}

