#' Confidence Intervals for Objects in tvReg
#'
#' confint is used to estimate the bootstrap confidence intervals for objects with class
#' attribute \code{tvlm}, \code{tvar}, \code{tvirf}, \code{tvsure}.
#' @param object Object of class \code{tvsure}, class \code{tvvar} or class \code{tvirf}.
#' @param parm A specification of which parameters are to be given confidence intervals, 
#' either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level Numeric, the confidence level required (between 0 and 1). 
#' @param runs (optional) Number of bootstrap replications.
#' @param tboot Type of wild bootstrap, choices 'wild'(default), 'wild2'. Option 'wild' uses the
#' distribution suggested by Mammen (1993) in the wild resampling, while 'wild2' uses the standard
#' normal.
#' @param ... Other parameters passed to specific methods.
#' @seealso \code{\link{tvLM}}, \code{\link{tvAR}}, \code{\link{tvVAR}},
#' \code{\link{tvSURE}}
#'
#' @return an object of class \code{tvsure} with BOOT, Lower and Upper different from NULL.
#'
#' @references
#' Chen, X. B., Gao, J., Li, D., and Silvapulle, P (2017) Nonparametric estimation and 
#' forecasting for time-varying coefficient realized volatility models,
#' \emph{Journal of Business \& Economic Statistics}, online, 1-13.
#'
#' Mammen, E (1993) Bootstrap and wild bootstrap for high dimensional linear models,
#' \emph{ Annals of Statistics}, 21, 255-285.
#'
#' @examples
#' \dontrun{
#' ##Calculation of confidence intervals for a TV-LM model
#' 
#' ##Generation of time-varying coefficients linear model
#' set.seed(42)
#' tau <- seq(1:200)/200
#' beta <- data.frame(beta1 = sin(2*pi*tau), beta2= 2*tau)
#' X1 <- rnorm(200)
#' X2 <- rchisq(200, df = 4)
#' error <- rt(200, df = 10)
#' y <- apply(cbind(X1, X2)*beta, 1, sum) + error
#' data <- data.frame(y = y, X1 = X1, X2 = X2)
#' 
#' ##Fitting the model and confidence interval calculation
#' model.tvlm <-  tvLM(y ~ 0 + X1 + X2, data = data, bw = 0.29)
#' tvci <- confint(model.tvlm, level = 0.95, runs = 20)
#' 
#' ##If a second confidence interval on the "same" object is calculated, 
#' ##for example with a different level, the calculation is faster
#' 
#' tvci.80 <- confint(tvci, level = 0.8)
#' }
#' @importFrom stats confint
#' @rdname confint.tvReg
#' @method confint tvlm
#' @export 
#'
confint.tvlm <- function(object, parm, level = 0.95, 
                         runs = 100, tboot = NULL, ...)
{
  if (!any(class(object) %in% c("tvlm", "tvar", "tvsure")))
    stop("\nConfidence intervals not implemented for this class.\n")
  if(runs <= 0)
    stop("\nVariable 'runs' accepts integer values greater than 0.\n")
  if (level <= 0 | level > 1)
    stop("\nVariable 'level' accepts values between 0 and 1.\n")
  if(is.null(tboot))
    tboot <- ifelse(!is.null(object$tboot), object$tboot, "wild")
  if(!(tboot %in% c("wild", "wild2")))
    tboot <- "wild"
  BOOT <- object$BOOT
  if(is.null(BOOT))
  {
    BOOT <- .tvboot(x = object, runs = runs, tboot = tboot)
  }
  else
  {
    if (object$tboot != tboot)
    {
      BOOT <- .tvboot(x = object, runs = runs, tboot = tboot)
    }
    else if (object$tboot == tboot & object$runs < runs)
    {
      temp <- .tvboot(x = object, runs = runs - object$runs, tboot = tboot)
      BOOT <- c(object$BOOT, temp)
    }
  }
  object$level <- level
  object$tboot <- tboot
  object$runs <- runs
  B <- object$tvcoef
  obs <- object$obs
  nvar <- ncol(object$x)
  if(any(class(object) == "tvsure"))
  {  
    neq <- object$neq
    nvar <- object$nvar
    obs <- object$obs
  }
  alpha <-  1 - level
  lower <- alpha/2
  upper <- 1 - lower
  mat.l <- matrix(0, nrow = obs, ncol = sum(nvar))
  mat.u <- matrix(0, nrow = obs, ncol = sum(nvar))
  temp <- matrix(NA, nrow = obs, ncol = runs)

  if(any(class(object) == "tvsure"))
  {
    for (m in 1:neq)
    {
      for (l in 1:nvar[m])
      {
        for (i in 1:runs)
        {
          temp[,i] <- BOOT[[i]][ , sum(nvar[1:m]) - nvar[m] + l]
        }
        sd.star <- apply(temp, 1, stats::sd)
        c.hat <- apply(abs(temp-B[ , sum(nvar[1:m]) - nvar[m] + l])/sd.star, 1, 
                       stats::quantile, prob = upper, na.rm = TRUE)
        mat.l[,sum(nvar[1:m]) - nvar[m] + l] <- B[, sum(nvar[1:m]) - nvar[m] + l] - c.hat*sd.star
        mat.u[,sum(nvar[1:m]) - nvar[m] + l] <- B[, sum(nvar[1:m]) - nvar[m] + l] + c.hat*sd.star
      }
    }
  }
  else
  {
    for (l in 1:nvar)
    {
      for (i in 1:runs)
      {
        temp[,i] <- BOOT[[i]][, l]
      }
      sd.star <- apply(temp, 1, stats::sd)
      c.hat <- apply(abs(temp-B[, l])/sd.star, 1, stats::quantile,
                     prob = upper, na.rm = TRUE)
      mat.l[, l] <- B[, l] - c.hat*sd.star
      mat.u[, l] <- B[, l] + c.hat*sd.star
    }
    colnames(mat.l) <- colnames(B)
    colnames(mat.u) <- colnames(B)
  }
  Lower <- mat.l
  Upper <- mat.u
  object$BOOT <- BOOT
  object$Lower <- Lower
  object$Upper <- Upper
  return(object)
}

#' @rdname confint.tvReg
#' @method confint tvar 
#' @export 
confint.tvar <- confint.tvlm

#' @rdname confint.tvReg
#' @method confint tvsure 
#' @export 
confint.tvsure <- confint.tvlm

#' @rdname confint.tvReg
#' @method confint tvirf
#' @export 
#'
confint.tvirf <- function(object, parm, level = 0.95, 
                          runs = 100, tboot = NULL , ...)
{
  if (class(object) != "tvirf")
    stop("\nConfidence intervals not implemented for this class.\n")
  if(runs <= 0)
    stop("\nVariable 'runs' accepts integer values greater than 0.\n")
  if (level <= 0 | level > 1)
    stop("\nVariable 'level' accepts values between 0 and 1.\n")
  if(is.null(tboot))
    tboot <- ifelse(!is.null(object$tboot), object$tboot, "wild")
  if(!(tboot %in% c("wild", "wild2")))
    tboot <- "wild"
  BOOT <- object$BOOT
  if(is.null(BOOT))
  {
    BOOT <- .tvboot(x = object, runs = runs, tboot = tboot)
  }
  else
  {
    if (object$tboot != tboot)
    {
      BOOT <- .tvboot(x = object, runs = runs, tboot = tboot)
    }
    else if (object$tboot == tboot & object$runs < runs)
    {
      temp <- .tvboot(x = object, runs = runs - object$runs, tboot = tboot)
      BOOT <- c(object$BOOT, temp)
    }
  }
  object$level <- level
  object$runs <- runs
  object$tboot <- tboot
  irf <- object$irf
  impulse <- object$impulse
  response <- object$response
  obs <- object$x$obs
  n.ahead <- object$n.ahead
  idx1 <- length(impulse)
  idx2 <- length(response)
  idx3 <- n.ahead + 1
  mat.l <- array(NA, dim = c(obs, idx2, n.ahead + 1))
  mat.u <- array(NA, dim = c(obs, idx2, n.ahead + 1))
  temp <- matrix(NA, nrow = obs, ncol = runs)
  alpha <-  1 - level
  lower <- alpha/2
  upper <- 1 - lower
  Lower <- list()
  Upper <- list()
  for (j in 1:idx1)
  {#impulse
    for (m in 1:idx2)
    {#response
      for (l in 1:idx3)
      {#horizon
        for (i in 1:runs)
        {
          temp[,i] <- BOOT[[i]][[j]][, m, l]
        }
        sd.star <- apply(temp, 1, stats::sd)
        if (sum(sd.star) == 0)
          c.hat <- apply(temp - irf[[j]][, m, l], 1, stats::quantile,
                         prob = upper, na.rm = TRUE)
        else
          c.hat <- apply(abs(temp - irf[[j]][, m, l])/sd.star, 1, stats::quantile,
                         prob = upper, na.rm = TRUE)
        mat.l[,m, l] <- irf[[j]][, m, l] - c.hat*sd.star
        mat.u[,m, l] <- irf[[j]][, m, l] + c.hat*sd.star
      }
    }
    dimnames(mat.l) <- list(NULL, response, NULL)
    dimnames(mat.u) <- list(NULL, response, NULL)
    Lower[[j]] <- mat.l
    Upper[[j]] <- mat.u
  }
  names(Lower) <- impulse
  names(Upper) <- impulse
  object$BOOT <- BOOT
  object$Lower <- Lower
  object$Upper <- Upper
  return(object)
}

