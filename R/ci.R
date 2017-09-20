
#' Confidence Intervals for Model Parameters of Objects in tvReg
#'
#' CI is used to estimate the bootstrap confidence intervals for objects with class
#' attribute \code{tvlm}, \code{tvar}, \code{tvirf}, \code{tvsure}.
#'
#' @param object Object of class \code{tvsure}, class \code{tvvar} or class \code{tvirf}.
#' @param level Numeric, the confidence level required (between 0 and 1). 
#' @param runs (optional) Number of bootstrap replications.
#' @param tboot Type of wild bootstrap, choices 'wild'(default), 'wild2'. Option 'wild' uses the
#' distribution suggested by Mammen (1993) in the wild resampling, while 'wild2' uses the standard
#' normal.
#' @param ... Other parameters passed to specific methods.
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
#' 
#' ## CI  for class 'tvlm'
#' tau <- seq(1:1000)/1000
#' beta <- data.frame(beta1 = sin(2*pi*tau), beta2= 2*tau)
#' X1 <- rnorm(1000)
#' X2 <- rchisq(1000, df = 4)
#' error <- rt(1000, df = 10)
#' y <- apply(cbind(X1, X2)*beta, 1, sum) + error
#' data <- data.frame(y = y, X1 = X1, X2 = X2)
#' model.tvlm <-  tvLM(y~0+X1+X2, data = data)
#' tvci <- CI(model.tvlm, level = 0.95, runs = 30)
#' tvci2 <- CI (tvci, level = 0.80)
#' tvci3 <- CI (tvci, level =0.80 , runs = 60)
#' plot(tvci)
#' 
#' ## CI for class 'tvsure'
#' data( "Kmenta", package="systemfit" )
#' eqDemand <- consump ~ price + income
#' eqSupply <- consump ~ price + farmPrice + trend
#' system <- list( demand = eqDemand, supply = eqSupply )
#'
#' tvfgls1.fit <- tvSURE(system, data = Kmenta, method="tvFGLS")
#'
#' ##Calculate 95% confidence interval of our object using a resampling of size 100
#' tvfgls1.fit <- CI (tvfgls1.fit, level = 0.95, runs = 100)
#'
#' ##Once the first confidence interval is calculated, the same resamples are used for other
#' ##confidence intervals, which makes the process faster for consecutives calls of function CI.
#' tvfgls1.fit <- CI (tvfgls1.fit, level = 0.90)
#' @export
CI <- function(object, ...) UseMethod("CI", object)


#' @rdname CI
#' @method CI default
#' @export
#'
CI.default <- function(object, level = 0, runs = 0, tboot = NULL, ...)
{
  if (class(object) != "tvlm" & class(object) != "tvar" & class(object) != "tvsure")
    stop("\nConfidence intervals not implemented for this class.\n")
  if (level < 0 | level > 1)
    stop("\nVariable 'level' accepts values between 0 and 1.\n")
  if (runs == 0 & object$runs != 0)
    runs <- object$runs
  if (runs == 0)
    runs <- ifelse(object$runs != 0, object$runs, 100)
  if(object$runs == 0)
    object$runs <- runs
  if(is.null(tboot))
    tboot <- ifelse(!is.null(object$tboot), object$tboot, "wild")
  if(!(tboot %in% c("wild", "wild2")))
    tboot <- "wild"
  if(is.null(object$tboot))
    object$tboot <- tboot
  if (level == 0)
    level <- ifelse(object$level != 0, object$level, 0.95)
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
  obs <- length(object$y)
  nvar <- ncol(object$x)
  if(class(object) == "tvsure")
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
  if(class(object) == "tvsure")
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
  }
  Lower <- mat.l
  Upper <- mat.u
  object$BOOT <- BOOT
  object$Lower <- Lower
  object$Upper <- Upper
  return(object)
}


#' @rdname CI
#' @method CI tvirf
#' @export
#'
CI.tvirf <- function(object, level = 0, runs = 0, tboot = NULL, ...)
{
  if (class(object) != "tvirf")
    stop("\nConfidence intervals not implemented for this class.\n")
  if (runs == 0)
    runs <- ifelse(object$runs != 0, object$runs, 100)
  if(object$runs == 0)
    object$runs <- runs
  if(is.null(tboot))
    tboot <- ifelse(!is.null(object$tboot), object$tboot, "wild")
  if(!(tboot %in% c("wild", "wild2")))
    tboot <- "wild"
  if(is.null(object$tboot))
    object$tboot <- tboot
  if (level == 0)
    level <- ifelse(object$level != 0, object$level, 0.95)
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
        mat.l[,m, l] <- irf[[j]][, m,l ] - c.hat*sd.star
        mat.u[,m, l] <- irf[[j]][, m,l ] + c.hat*sd.star
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

