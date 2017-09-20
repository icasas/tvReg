##' @rdname tvReg-internals
##' @keywords internal
.tvIRF <- function (x, impulse, response, y.names, n.ahead, ortho, cumulative,
                    ortho.cov, bw.cov)
{
  if (class(x) == "tvvar")
  {
    if (ortho)
    {
      result <-  tvPsi(x, nstep = n.ahead, ortho.cov, bw.cov = bw.cov)
      irf <- result$Psi
      bw.cov <- result$bw.cov
    }
    else
    {
      irf <- tvPhi(x, nstep = n.ahead)
    }
  }
  dimnames(irf) <- list(NULL, y.names, y.names, NULL)
  idx <- length(impulse)
  irs <- list()
  for (i in 1:idx)
  {
    irs[[i]] <- irf[, response, impulse[i], 1:(n.ahead + 1)]
    #i have to figure out the names here
    if (cumulative)
    {
      if (length(response) > 1)
        irs[[i]] <- apply(irs[[i]], 2, cumsum)
      if (length(response) == 1)
      {
        tmp <- matrix(cumsum(irs[[1]]))
        colnames(tmp) <- response
        irs[[1]] <- tmp
      }
    }
  }
  names(irs) <- impulse
  result <- irs
  if(ortho)
    return(list(irf = result, bw.cov = bw.cov))
  return(list(irf = result, bw.cov = NULL))
}


##' @rdname tvReg-internals
##' @keywords internal
.tvOLS.cv <- function(bw, x, y, z = NULL, est = "lc", tkernel = "Epa", singular.ok = TRUE)
{
  x <- as.matrix(x)
  obs <- nrow(x)
  fitted <- resid.2 <- numeric(obs)
  nvar <- ncol(x)
  bw <- abs(bw)
  if(!is.null(z))
    grid <- z
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
    result <- stats::lm.wfit(x = xtemp, y = ytemp, w = kernel.bw[myindex],
                             singular.ok = singular.ok)
    theta[t,] <- result$coef[1:nvar]
    fitted[t] <- crossprod(x[t, !is.na(theta[t,])], theta[t, !is.na(theta[t,])])
  }
  resid.2 <- (y - fitted)^2
  return(mean(resid.2))
}

##' Create List of Control Parameters for tvSURE
##'
##' Create a list of control pararameters for function \code{\link{tvSURE}}.
##' All control parameters that are not passed to this function are set to
##' default values.
##'
##' If the estimation is iterative FGLS with \code{maxiter}>1, the convergence criterion is
##' \deqn{\sqrt{ \frac{ \sum_{i, j}
##' (B_{i,j,g} - B{i, j, g-1})^2 }{ \sum_{i, j} B_{i, j, g-1}^2 }} < \code{tol}}
##' (\eqn{B_{i, j,g}} is the ith, jth coefficient of the gth iteration step).
##' @rdname tvsystem-internals
##' @param maxiter maximum number of iterations for the iterative FGLS
##' estimations.
##' @param tol tolerance level indicating when to stop the iteration for the iterative FGLS estimations
##' @return A list of the above components.
##'
##' @export
tvsure.control <- function(maxiter = 1, tol = 1e-05)
{
  result <- list()
  if (maxiter <= 0 || round(maxiter) != maxiter) {
    stop("control parameter 'maxiter' must be a positive integer\n")
  }
  result$maxiter <- maxiter
  if (tol <= 0 || !is.numeric(tol) || length(tol) != 1) {
    stop("control parameter 'tol' must be a positive scalar\n")
  }
  result$tol <- tol
  return(result)
}

