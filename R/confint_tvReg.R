#' Confidence Intervals for Objects in tvReg
#'
#' confint is used to estimate the bootstrap confidence intervals for objects with class
#' attribute \code{tvlm}, \code{tvar}, \code{tvirf}, \code{tvsure} and \code{tvplm}.
#' @param object An object used to select a method.
#' @param parm A specification of which parameters are to be given confidence intervals, 
#' either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level Numeric, the confidence level required (between 0 and 1). 
#' @param runs (optional) Number of bootstrap replications.
#' @param tboot Type of wild bootstrap, choices 'wild'(default), 'wild2' and 'lbwb'. 
#' Option 'wild' uses the distribution suggested by Mammen (1993) in the wild resampling, 
#' while 'wild2' uses the standard normal. The option of Local Blockwise Wild Bootstrap, 'lbwb', 
#' uses block-bootstrap suggested by Lin et al (2025).
#' @param ... Other parameters passed to specific methods.
#' @seealso \code{\link{tvLM}}, \code{\link{tvAR}}, \code{\link{tvVAR}},
#' \code{\link{tvSURE}}
#'
#' @return an object of class \code{tvsure} with BOOT, Lower and Upper different from NULL.
#'
#' @references
#' Chen, X. B., Gao, J., Li, D., and Silvapulle, P. (2017) Nonparametric estimation and 
#' forecasting for time-varying coefficient realized volatility models,
#' \emph{Journal of Business and Economic Statistics}, 36, 88-100.
#'
#' Mammen, E. (1993) Bootstrap and wild bootstrap for high dimensional linear models,
#' \emph{ Annals of Statistics}, 21, 255-285.
#' 
#' Lin, Y., Song, M., and van der Sluis, B. (2025) Bootstrap Inference for Linear 
#' Time-Varying Coefficient Models in Locally Stationary Time Series, 
#' \emph{Journal of Computational and Graphical Statistics}, 34:2, 654-667.
#'
#' @examples
#' \dontrun{
#' ##Calculation of confidence intervals for a TVLM model
#' 
#' ##Generation of time-varying coefficients linear model
#' set.seed(42)
#' tau <- seq(1:200)/200
#' beta <- data.frame(beta1 = sin(2*pi*tau), beta2= 2*tau)
#' X1 <- stats::rnorm(200)
#' X2 <- rchisq(200, df = 4)
#' error <- rt(200, df = 10)
#' y <- apply(cbind(X1, X2)*beta, 1, sum) + error
#' data <- data.frame(y = y, X1 = X1, X2 = X2)
#' 
#' ##Fitting the model and confidence interval calculation
#' model.tvlm <-  tvLM(y ~ 0 + X1 + X2, data = data, bw = 0.29)
#' tvci <- confint(model.tvlm, level = 0.95, runs = 20, tboot = "lbwb")
#' 
#' ##If a second confidence interval on the "same" object is calculated, 
#' ##for example with a different level, the calculation is faster
#' tvci.80 <- confint(tvci, level = 0.8)
#' }
#' @importFrom stats confint
#' @rdname confint.tvReg
#' @method confint tvlm
#' @export 
#'


confint.tvlm <- function(object, parm, level = 0.95, 
                         runs = 100, tboot = c("wild", "wild2", "lbwb"), ...)
{
  if (!any(class(object) %in% c("tvlm", "tvar", "tvsure")))
    stop("\nConfidence intervals not implemented for this class.\n")
  if(runs <= 0)
    stop("\nVariable 'runs' accepts integer values greater than 0.\n")
  if (level <= 0 | level > 1)
    stop("\nVariable 'level' accepts values between 0 and 1.\n")
  tboot <- match.arg(tboot)
  BOOT <- object$BOOT
  if(is.null(BOOT))
  {
    boot_out <- .tvboot(x = object, runs = runs, tboot = tboot)
    BOOT <- boot_out$BOOT
    object$beta.pilot <- boot_out$beta.pilot
  }
  else
  {
    if (object$tboot != tboot)
    {
      boot_out <- .tvboot(x = object, runs = runs, tboot = tboot)
      BOOT <- boot_out$BOOT
      object$beta.pilot <- boot_out$beta.pilot
    }
    else if (object$tboot == tboot & object$runs < runs)
    {
      boot_out <- .tvboot(x = object, runs = runs - object$runs, tboot = tboot)
      BOOT <- c(object$BOOT, boot_out$BOOT)
      object$beta.pilot <- boot_out$beta.pilot
    }
  }
  object$level <- level
  object$tboot <- tboot
  object$runs <- runs
  B <- object$coefficients
  if(tboot == "lbwb")
  {
    B.pilot <- object$beta.pilot
    if(is.null(B.pilot)) stop("Pilot coefficients are missing for LBWB.")
  }
  obs <- object$obs
  nvar <- ncol(object$x)
  if(inherits(object, "tvsure"))
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
  if(inherits(object, "tvsure"))
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
        if (tboot == "lbwb") 
        {
          B.ref <- B.pilot[, sum(nvar[1:m]) - nvar[m] + l]
        } 
        else 
        {
          B.ref <- B[, sum(nvar[1:m]) - nvar[m] + l]
        }
        
        c.hat <- apply(abs(temp - B.ref)/sd.star, 1, 
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
      # Replace your current B.ref line with this:
      if (tboot == "lbwb") {
        B.ref <- B.pilot[, l]
      } else {
        B.ref <- B[, l]
      }
      
      c.hat <- apply(abs(temp - B.ref)/sd.star, 1, 
                     stats::quantile, prob = upper, na.rm = TRUE)
      mat.l[, l] <- B[, l] - c.hat*sd.star
      mat.u[, l] <- B[, l] + c.hat*sd.star
    }
    colnames(mat.l) <- colnames(B)
    colnames(mat.u) <- colnames(B)
  }
  object$BOOT <- BOOT
  object$Lower <- mat.l
  object$Upper <- mat.u
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
#' @method confint tvvar 
#' @export 
confint.tvvar <- function(object, parm, level = 0.95, 
                          runs = 100, tboot = c("wild", "wild2", "lbwb"),
                          type = c("pointwise", "simultaneous"), ...)
{
  stop("\nConfidence intervals are calculated on the 'tvirf' object, not directly on 'tvvar'. Use confint(tvirf_object) instead.\n")
}


#' @rdname confint.tvReg
#' @method confint tvirf
#' @export 
#'
confint.tvirf <- function(object, parm, level = 0.95, 
                          runs = 100, tboot = c("wild", "wild2", "lbwb"),
                          type = c("pointwise", "simultaneous"), ...)
{
  if (!inherits(object, "tvirf"))
    stop("\nConfidence intervals not implemented for this class.\n")
  if (runs <= 0)
    stop("\nVariable 'runs' accepts integer values greater than 0.\n")
  if (level <= 0 | level > 1)
    stop("\nVariable 'level' accepts values between 0 and 1.\n")
  
  tboot <- match.arg(tboot)
  type  <- match.arg(type)
  
  BOOT <- object$BOOT
  if (is.null(BOOT) || is.null(object$tboot) || (object$tboot != tboot)) 
  {
    boot_out <- .tvboot(x = object, runs = runs, tboot = tboot)
    BOOT <- boot_out$BOOT
    object$beta.pilot <- boot_out$beta.pilot
  } 
  else if (object$runs < runs) 
  {
    boot_out <- .tvboot(x = object, runs = runs - object$runs, tboot = tboot)
    BOOT <- append(object$BOOT, boot_out$BOOT)
    object$beta.pilot <- boot_out$beta.pilot
  }
  
  object$level <- level
  object$runs  <- runs
  object$tboot <- tboot
  object$type  <- type
  
  irf <- object$irf
  if (tboot == "lbwb") 
  {
    if (is.null(object$beta.pilot)) stop("Pilot coefficients missing for LBWB.")
    if (inherits(object$x, "tvvar")) 
    {  
      VAR_pilot <- object$x
    }
    else 
    {
      VAR_pilot <- eval.parent(object$x)
    }
    VAR_pilot$coefficients <- object$beta.pilot 
    
    irf_ref <- .tvIRF(x = VAR_pilot, impulse = object$impulse, 
                      response = object$response, y.names = colnames(VAR_pilot$y),
                      n.ahead = object$n.ahead, ortho = object$ortho, 
                      cumulative = object$cumulative, 
                      ortho.cov = object$ortho.cov, bw.cov = object$bw.cov, 
                      unit.shock = object$unit.shock)$irf
  } 
  else
  {
    irf_ref <- irf
  }
  
  idx1 <- length(object$impulse)
  idx2 <- length(object$response)
  idx3 <- object$n.ahead + 1
  alpha <- 1 - level
  upper_prob <- 1 - (alpha / 2)
  Lower <- Upper <- list()
  
  nobs <- object$x$obs
  for (j in 1:idx1) 
  {
    mat.l <- array(NA, dim = c(nobs, idx2, idx3))
    mat.u <- array(NA, dim = c(nobs, idx2, idx3))
    
    for (m in 1:idx2)
    {
      for (l in 1:idx3) 
      {
        temp <- sapply(BOOT, function(b) b[[j]][, m, l])
        if (is.null(dim(temp)))
          temp <- matrix(temp, nrow = nobs, ncol = length(BOOT))
        
        sd.star <- apply(temp, 1, stats::sd)
        sd.star[sd.star == 0] <- 1e-8 # Avoid division by zero
        center <- irf_ref[[j]][, m, l]
        
        t.stat <- abs(temp - center) / sd.star
        
        if (type == "pointwise") 
        {
          c.hat <- apply(t.stat, 1, stats::quantile, prob = upper_prob, na.rm = TRUE)
        } 
        else 
        {
          sup_t <- apply(t.stat, 2, max, na.rm = TRUE)
          c.hat <- stats::quantile(sup_t, prob = level, na.rm = TRUE)
        }
        
        mat.l[, m, l] <- irf[[j]][, m, l] - c.hat * sd.star
        mat.u[, m, l] <- irf[[j]][, m, l] + c.hat * sd.star
      }
    }
    dimnames(mat.l) <- dimnames(mat.u) <- list(NULL, object$response, NULL)
    Lower[[j]] <- mat.l
    Upper[[j]] <- mat.u
  }
  
  names(Lower) <- names(Upper) <- object$impulse
  object$BOOT  <- BOOT
  object$Lower <- Lower
  object$Upper <- Upper
  return(object)
}


#' @rdname confint.tvReg
#' @method confint tvplm
#' @export 
#'
confint.tvplm <- function(object, parm, level = 0.95,  
                          runs = 100, tboot = c("wild", "wild2", "lbwb"),
                          type = c("pointwise", "simultaneous"), ...)
{
  if (!inherits(object, "tvplm"))
    stop("\nConfidence intervals not implemented for this class.\n")
  if (runs <= 0)
    stop("\nVariable 'runs' accepts integer values greater than 0.\n")
  if (level <= 0 | level > 1)
    stop("\nVariable 'level' accepts values between 0 and 1.\n")
  
  tboot <- match.arg(tboot)
  type  <- match.arg(type)
  
  BOOT <- object$BOOT
  if (is.null(BOOT) || is.null(object$tboot) || (object$tboot != tboot))
  {
    boot_out <- .tvboot(x = object, runs = runs, tboot = tboot)
    BOOT <- boot_out$BOOT
    object$beta.pilot <- boot_out$beta.pilot
  }
  else if (object$runs < runs)
  {
    boot_out <- .tvboot(x = object, runs = runs - object$runs, tboot = tboot)
    BOOT <- c(object$BOOT, boot_out$BOOT)
    object$beta.pilot <- boot_out$beta.pilot
  }
  
  object$level <- level
  object$runs  <- runs
  object$tboot <- tboot
  object$type  <- type
  
  obs  <- object$obs
  nvar <- object$nvar
  B    <- object$coefficients
  
  if (tboot == "lbwb")
  {
    B.pilot <- object$beta.pilot
    if (is.null(B.pilot)) stop("Pilot coefficients missing for LBWB.")
  }
  
  mat.l <- matrix(NA, nrow = obs, ncol = nvar)
  mat.u <- matrix(NA, nrow = obs, ncol = nvar)
  temp  <- matrix(NA, nrow = obs, ncol = runs)
  
  alpha <- 1 - level
  upper <- 1 - (alpha / 2)
  
  for (l in 1:nvar) 
  {
    for (i in 1:runs) 
    {
      temp[, i] <- BOOT[[i]][, l]
    }
    sd.star <- apply(temp, 1, stats::sd)
    sd.star[sd.star == 0] <- 1e-8
    
    B.ref <- if (tboot == "lbwb") B.pilot[, l] else B[, l]
    t.stat <- abs(temp - B.ref) / sd.star
    
    if (type == "pointwise") 
    {
      c.hat <- apply(t.stat, 1, stats::quantile, prob = upper, na.rm = TRUE)
    }
    else 
    {
      sup_t <- apply(t.stat, 2, max, na.rm = TRUE)
      c.hat <- stats::quantile(sup_t, prob = level, na.rm = TRUE)
    }
    
    mat.l[, l] <- B[, l] - c.hat * sd.star
    mat.u[, l] <- B[, l] + c.hat * sd.star
  }
  
  colnames(mat.l) <- colnames(B)
  colnames(mat.u) <- colnames(B)
  
  object$BOOT  <- BOOT
  object$Lower <- mat.l
  object$Upper <- mat.u
  return(object)
}


#' @rdname tvReg-internals
#' @keywords internal
.tvLM.ci <- function(x, yboot)
{
  obs <- NROW(yboot)
  runs <- NCOL(yboot)
  bw <- x$bw
  z <- x$z
  ez <- x$ez
  if(!is.null(z))
    grid <- z
  else
    grid <- (1:obs)/obs
  if (is.null(ez))
    ez <- grid
  tkernel <- x$tkernel
  est <- x$est
  nvar <- NCOL(x$x)
  eobs <- NROW(ez)
  theta <- matrix(0, eobs, nvar)
  x <- x$x
  BOOT <- vector("list", runs)
  for (t in 1:eobs)
  { 
    tau0 <- grid - ez[t]
    kernel.bw <- .kernel(x = tau0, bw = bw, tkernel = tkernel)
    k.index <- which(kernel.bw != 0)
    xtemp <- x[k.index, ]
    if (est == "ll")
      xtemp <- cbind(xtemp, xtemp * tau0[k.index])
    for (k in 1:runs)
    {
      result <- stats::lm.wfit(x = as.matrix(xtemp), y = yboot[k.index, k], w = kernel.bw[k.index])
      BOOT[[k]] <-rbind(BOOT[[k]], result$coefficients[1:nvar])
    }
  }
  return(BOOT)
}


#' @rdname tvReg-internals
#' @keywords internal
.tvSURE.ci <- function(x, yboot)
{
  obs <- NROW(yboot)
  runs <- length(yboot)
  method <- x$method
  tkernel <- x$tkernel
  obs <- x$obs
  neq <- x$neq
  Sigma <- x$Sigma
  bw <- x$bw
  bw.cov <- x$bw.cov
  z <- x$z
  ez <- x$ez
  R <- x$R
  r <- x$r
  if(!is.null(z))
    grid <- z
  else
    grid <- (1:obs)/obs
  ez <- grid
  tkernel <- x$tkernel
  est <- x$est
  nvar <- x$nvar
  eobs <- NROW(ez)
  theta <- matrix(0, eobs, sum(nvar))
  x <- x$x
  BOOT = resid <- vector("list", runs)
  if (length(bw) == 1)
    bw <- rep(bw, neq)
  if(!is.null(R))
  {
    R <- as.matrix(R)
    if (is.null(r))
      r <- rep(0, NROW(R))
    else if (length(r) == 1)
      r <- rep(r, NROW(R))
  }
  if(method %in% c("identity", "tvOLS", "tvFGLS"))
  {
    Sigma <- array(rep(diag(1, neq), obs), dim = c(neq, neq, obs))
  }
  for (t in 1:eobs)
  {
    tau0 <- grid - ez[t]
    y.kernel <- NULL
    x.kernel <- vector ("list", neq)
    eSigma <- eigen(Sigma[,,t], TRUE)
    if (any(eSigma$value <= 0))
      stop("\n'Sigma' is not positive definite.\n")
    A <- diag(eSigma$values^-0.5) %*% t(eSigma$vectors) %x% Matrix::Diagonal(obs)
    for (i in 1:neq)
    {
      mykernel <- sqrt(.kernel(x = tau0, bw = bw[i], tkernel = tkernel))          
      xtemp <- x[[i]] * mykernel
      if(est == "ll")
        xtemp <- cbind(xtemp, xtemp * tau0)
      x.kernel[[i]] <- xtemp
    }
    x.star <- A %*% Matrix::bdiag(x.kernel)
    s0 <- Matrix::crossprod(x.star)
    if(!is.null(R))
    {
      Sinv <- qr.solve(s0)
      mat <- Sinv %*% t(R) %*% qr.solve(R %*% Sinv %*% t(R))
    }
    for(k in 1:runs)
    {
      y.kernel <- NULL
      for(i in 1:neq)
      {
        mykernel <- sqrt(.kernel(x = tau0, bw = bw[i], tkernel = tkernel))
        y.kernel <- cbind(y.kernel, yboot[[k]][, i] * mykernel)
      }
      y.star <- A %*% methods::as(y.kernel, "sparseVector")
      T0 <- Matrix::crossprod(x.star, y.star)
      result <- try(qr.solve(s0, T0), silent = TRUE)
      if(inherits(result, "try-error"))
        result <- try(qr.solve(s0, T0, tol = .Machine$double.xmin ), silent = TRUE)
      if(inherits(result, "try-error"))
        stop("\nSystem is computationally singular, the inverse cannot be calculated. 
             Possibly, the 'bw' is too small for values in 'ez'.\n")
      if (est == "ll") 
      {
        temp <- numeric(sum(nvar))
        start.ind <- 1
        start.res <- 1
        for (i in 1:neq) 
        {
          temp[start.ind:(start.ind + nvar[i] - 1)] <- result[start.res:(start.res + nvar[i] - 1)]
          start.ind <- start.ind + nvar[i]
          start.res <- start.res + 2 * nvar[i]
        }
        theta[t, ] <- temp
      } 
      else 
      {
        theta[t, ] <- result
      }
      if(!is.null(R))
      {
        theta[t, ] <- drop(theta[t, ] - mat %*% (R %*% theta[t, ] - r))
      }
      BOOT[[k]] <- rbind (BOOT[[k]], theta[t,])
      resid[[k]] <- rbind(resid[[k]], yboot[[k]][t,] - as.numeric(Matrix::bdiag(lapply(x,"[",t, ,drop = FALSE))%*%theta[t,]))
    }
  }
  if(method == "tvFGLS")
  {
    BOOT = Cov <- vector("list", runs)
    for (k in 1:runs)
    {
      Cov[[k]] <- tvCov(x = resid[[k]], bw = bw.cov, tkernel = tkernel)
    }
    for (t in 1:eobs)
    {
      tau0 <- grid - ez[t]
      y.kernel <- NULL
      x.kernel <- vector ("list", neq)
      for (i in 1:neq)
      {
        mykernel <- sqrt(.kernel(x = tau0, bw = bw[i], tkernel = tkernel))          
        xtemp <- x[[i]] * mykernel
        if(est == "ll")
          xtemp <- cbind(xtemp, xtemp * tau0)
        x.kernel[[i]] <- xtemp
      }
      for(k in 1:runs)
      {
        eSigma <- eigen(Cov[[k]][,,t], TRUE)
        if (any(eSigma$value <= 0))
        {
          warning("\n'Sigma' is not positive definite, re-run fitting.\n")
          bw.cov2 <- bwCov(x = resid[[k]], cv.block = floor(obs/10), tkernel = tkernel)
          Cov[[k]] <- tvCov(x = resid[[k]], bw = bw.cov2 , tkernel = tkernel)
          eSigma <- eigen(Cov[[k]][,,t], TRUE)
        }
        A <- diag(eSigma$values^-0.5) %*% t(eSigma$vectors) %x% Matrix::Diagonal(obs)
        x.star <- A %*% Matrix::bdiag(x.kernel)
        s0 <- Matrix::crossprod(x.star)
        y.kernel <- NULL
        for(i in 1:neq)
        {
          mykernel <- sqrt(.kernel(x = tau0, bw = bw[i], tkernel = tkernel))
          y.kernel <- cbind(y.kernel, yboot[[k]][, i] * mykernel)
        }
        y.star <- A %*% methods::as(y.kernel, "sparseVector")
        T0 <- Matrix::crossprod(x.star, y.star)
        result <- try(qr.solve(s0, T0), silent = TRUE)
        if(inherits(result, "try-error"))
          result <- try(qr.solve(s0, T0, tol = .Machine$double.xmin ), silent = TRUE)
        if(inherits(result, "try-error"))
          stop("\nSystem is computationally singular, the inverse cannot be calculated. 
               Possibly, the 'bw' is too small for values in 'ez'.\n")
        if (est == "ll") 
        {
          temp <- numeric(sum(nvar))
          start.ind <- 1
          start.res <- 1
          for (i in 1:neq) 
          {
            temp[start.ind:(start.ind + nvar[i] - 1)] <- result[start.res:(start.res + nvar[i] - 1)]
            start.ind <- start.ind + nvar[i]
            start.res <- start.res + 2 * nvar[i]
          }
          theta[t, ] <- temp
        } 
        else 
        {
          theta[t, ] <- result
        }
        if(!is.null(R))
        {
          Sinv <- qr.solve(s0)
          theta[t, ] <- drop(theta[t, ] - Sinv %*% t(R) %*% qr.solve(R %*% Sinv %*% t(R)) %*% (R %*% theta[t, ] - r))
        }
        BOOT[[k]] <- rbind (BOOT[[k]], theta[t,])
      }
    }
  }
  return(BOOT)
}

.tvPLM.ci<-function(x, yboot)
{
  obs <- x$obs
  runs <- NCOL(yboot)
  bw <- x$bw
  z <- x$z
  ez <- x$ez
  if(!is.null(z))
    grid <- z
  else
    grid <- (1:obs)/obs
  if (is.null(ez))
    ez <- grid
  tkernel <- x$tkernel
  est <- x$est
  neq <- x$neq
  nvar <- x$nvar
  eobs <- NROW(ez)
  method <- x$method
  xnew <- matrix(NA, neq*eobs, nvar)
  if(x$est == "ll")
    xnew <- matrix(NA, neq*eobs, 2*nvar)
  ynew <- yboot[, 1]
  BOOT <- vector("list", runs)
  D <- t(cbind(rep(-1, neq-1), diag(1, neq-1)))%x%rep(1, obs)
  for (t in 1:eobs)
  { 
    tau0 <- grid - ez[t]
    xtemp <- x$x
    if (est=="ll")
      xtemp <- cbind(xtemp, xtemp * tau0)
    if (method != "within")
    {
      kernel.bw <- sqrt(.kernel(tau0, bw, tkernel))
      k.index <- which(kernel.bw != 0)
      if (length (k.index) < 3)
        stop("Bandwidth is too small for values in 'ez'.\n")
      if(method == "pooling")
        Sigma <- diag(1, eobs)
      else
        Sigma <- x$Sigma
      eSigma <- eigen(Sigma, TRUE)
      if (any(eSigma$value <= 0))
        stop("\n'Sigma' is not positive definite.\n")
      A <- diag(eSigma$values^-0.5) %*% t(eSigma$vectors)
      for (i in 1:neq)
      {
        ind <- (i-1)*eobs + (1:eobs)
        xnew[ind, ] <- A%*%(xtemp[ind,] * kernel.bw)
      }
      mat.inv <- qr.solve(crossprod(xnew))
    }
    else
    {
      kernel.bw <- .kernel(tau0, bw, tkernel)
      k.index <- which(kernel.bw != 0)
      if (length (k.index) < 3)
        stop("Bandwidth is too small for values in 'ez'.\n")
      WH <- diag(neq)%x%diag(kernel.bw)
      x.tilde <- crossprod(xtemp, WH)
      temp <- Matrix::crossprod(D, WH)
      MH <- diag(neq*obs) - D %*% qr.solve(temp%*%D)%*%temp
      SH <- Matrix::crossprod(MH, WH)%*%MH
      x.tilde <- Matrix::crossprod(xtemp, SH)
      s0 <- x.tilde %*% xtemp
    }
    for (k in 1:runs)
    {
      if(method != "within")
      {
        for (i in 1:neq)
        {
          ind <- (i-1)*eobs + (1:eobs)
          ynew[ind] <- A%*%(yboot[ind, k] * kernel.bw)
        }
        result <- mat.inv%*%crossprod(xnew, ynew)
      }
      else
      {
        T0 <- x.tilde %*% yboot[, k]
        result <- try(qr.solve(s0, T0), silent = TRUE)
        if(inherits(result, "try-error"))
          result <- try(qr.solve(s0, T0, tol = .Machine$double.xmin), silent = TRUE)
      }
      BOOT[[k]] <-rbind(BOOT[[k]], result[1:nvar])
    }
  }
  return(BOOT)
}



