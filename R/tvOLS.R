#' Time-Varying Ordinary Least Squares
#'
#' \code{tvOLS} is used to fit univariate linear models with time-varying coefficients.
#'
#' @param x A matrix with all regressors.
#' @param y A vector with dependent variable.
#' @param z A vector with the variable over which coefficients are smooth over.
#' @param bw A numeric vector.
#' @inheritParams tvSURE
#' @param singular.ok	Logical. If FALSE, a singular model is an error.
#' @return A list with the estimates, fitted and residuals values.
#'
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
#' @seealso \code{\link{bw}} for bandwidth selection and \code{\link{tvLM}}
#' @export

tvOLS <- function(x, y, z = NULL, bw, est = c("lc", "ll"), 
                  tkernel = c("Epa", "Gaussian"), singular.ok = singular.ok)
{
  x <- as.matrix(x)
  y <- as.numeric(y)
  obs <- NROW(x)
  if(!identical(length(y), obs))
    stop("\nDimensions of 'x' and 'y' are not compatible.\n")
  if(!is.numeric(bw))
    stop ("Parameter 'bw' should be a scalar. \n")
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  if(!(tkernel %in% c("Epa", "Gaussian")))
    tkernel <- "Epa"
  if(!(est %in% c("lc", "ll")))
    est <- "lc"
  fitted <- numeric(obs)
  resid <- numeric(obs)
  nvar <- NCOL(x)
  theta <- matrix(0, obs, nvar)
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
    x2 <- grid - grid[t]
    kernel.bw <- .kernel(x = x2, bw = bw, tkernel = tkernel, N = 1)
    myindex <- which(kernel.bw != 0)
    if (length (myindex) < 3)
      stop("\nBandwidth is too small.\n.")
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
  resid <- y - fitted
  return(list(tvcoef = theta, fitted = fitted, residuals = resid))
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
