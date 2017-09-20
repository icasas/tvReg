##' Time-Varying Ordinary Least Squares
##'
##' \code{tvOLS} is used to fit univariate linear models with time-varying coefficients.
##'
##' @param x A matrix with all regressors.
##' @param y A vector with dependent variable.
##' @param z A vector with the variable over which coefficients are smooth over.
##' @param bw A numeric vector.
##' @inheritParams tvSURE
##' @param singular.ok	Logical. If FALSE, a singular model is an error.
##' @return A list with the estimates, fitted and residuals values.
##'
##' @examples
##' tau <- seq(1:1000)/1000
##' beta <- data.frame(beta1 = sin(2*pi*tau), beta2= 2*tau)
##' X <- data.frame(X1= rnorm(1000), X2 = rchisq(1000, df = 4))
##' error <- rt(1000, df = 10)
##' y <- apply(X*beta, 1, sum) + error
##' coef.lm <- stats::lm(y~0+X1+X2, data = X)$coef
##' bandw <- bw (x = X, y = y)
##' coef.tvlm <-  tvOLS(x = X, y = y, bw = bandw)$tvcoef
##' plot(tau,beta[, 1], type="l", main="", ylab = expression(beta[1]), xlab = expression(tau),
##' ylim = range(beta[,1], coef.tvlm[, 1]))
##' abline(h = coef.lm[1], col = 2)
##' lines(tau, coef.tvlm[, 1], col = 4)
##' legend("topright", c(expression(beta[1]), "lm", "tvlm"), col = c(1, 2, 4), bty="n", lty = 1)
##'
##' @seealso \code{\link{bw}} for bandwidth selection and \code{\link{tvLM}}
##' @export

tvOLS <- function(x, y, z = NULL, bw, est = "lc", tkernel = "Epa",
                  singular.ok = singular.ok)
{
  x <- as.matrix(x)
  y <- as.numeric(y)
  if(!is.numeric(bw))
    stop ("Parameter bw should be a scalar. \n")
  if(tkernel != "Epa" & tkernel != "Gaussian")
    tkernel <- "Epa"
  if(est != "lc" & est != "ll")
    est <- "lc"
  obs <- nrow(x)
  fitted <- numeric(obs)
  resid <- numeric(obs)
  nvar <- ncol(x)
  theta <- matrix(0, obs, nvar)
  if(!is.null(z))
    grid <- z
  else
    grid <- (1:obs)/obs
  for (t in 1:obs)
  {
    x2 <- grid - grid[t]
    kernel.bw <- .kernel(x = x2, bw = bw, tkernel = tkernel, N = 1)
    myindex <- which(kernel.bw != 0)
    if (sum(myindex) == 0)
      return (.Machine$double.xmax)
    xtemp <- x[myindex, ]
    if (est=="ll")
      xtemp <- cbind(xtemp, xtemp * x2[myindex])
    ytemp <- y[myindex]
    result <- stats::lm.wfit(x = xtemp, y = ytemp, w = kernel.bw[myindex])
    theta[t,] <- result$coef[1:nvar]
    fitted[t] <- crossprod(x[t, !is.na(theta[t,])], theta[t, !is.na(theta[t,])])
  }
  resid <- y - fitted
  return(list(tvcoef = theta, fitted = fitted, residuals = resid))
}
