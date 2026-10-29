#' Bandwidth Selection by Cross-Validation
#'
#' Calculate bandwidth(s) by cross-validation for functions \code{tvSURE}, \code{tvVAR} and \code{tvLM}.
#'
#' @param x An object used to select a method.
#' @param ... Other parameters passed to specific methods.
#' @param y A matrix or vector with the dependent variable(s).
#' @param z A vector with the variable over which coefficients are smooth over.
#' @param cv.block A positive scalar with the size of the block in leave-one block-out cross-validation.
#' By default 'cv.block=0' meaning leave-one-out cross-validation.
#' @param lmcv Logical. If \code{FALSE} (default), general cross-validation is used. 
#' If \code{TRUE}, local modified cross-validation (LMCV) is applied.
#' @param local.bw A numeric value or \code{NULL}. The window size used to weight 
#' the local modified cross-validation. If \code{NULL}, it is automatically calculated.
#' @param est The nonparametric estimation method, one of "lc" (default) for linear constant
#' or "ll" for local linear.
#' @param tkernel A character, either "Triweight" (default), "Epa" or "Gaussian" kernel function.
#' @param singular.ok Logical. If FALSE, a singular model is an error.
#' @param method A character with the choice of panel model/estimation method:
#' \itemize{
#'   \item \code{"tvPOLS"} (default): the data is pooled and estimated with time-varying OLS. 
#'   No individual or time effects are estimated.
#'   \item \code{"tvFE"}: individual effects which might be correlated with 
#'   the regressors are estimated.
#'   \item \code{"tvRE"}: individual effects are considered random and independent
#'   of the regressors.
#' }
#'
#' @return A scalar or a vector of scalars with the bandwidth(s) selected by cross-validation.
#'
#' @examples
#' ##Generate data
#' tau <- seq(1:200)/200
#' beta <- data.frame(beta1 = sin(2*pi*tau), beta2 =  2*tau)
#' X <- data.frame(X1 = stats::rnorm(200), X2 = rchisq(200, df = 4))
#' error <- rt(200, df = 10)
#' y <- apply(X*beta, 1, sum) + error
#' 
#' ##Select bandwidth by cross-validation
#' bw <- bw(X, y, est = "ll", tkernel = "Gaussian")
#'
#' @export
bw <- function(x, ...)  UseMethod("bw", x)

#' @rdname bw
#' @method bw default
#' @export
bw.default <- function(x, y, z = NULL, cv.block = 0, lmcv = FALSE,
                       local.bw = NULL, est = c("lc", "ll"),  
                       tkernel = c("Triweight", "Epa", "Gaussian"),
                       singular.ok = TRUE, ...)
{ 
  if(!inherits(x, c("matrix", "data.frame", "vector", "numeric", "integer")))
    stop("'x' should be a matrix, a vector or a data frame. \n")
  if(is.null(y))
    stop("Parameter 'y' missing. \n")
  if(sum(is.na(y))>0 | sum(is.na(x))>0)
    stop("There are NA values in your data, please enter only complete cases. \n")
  if (is.logical(local.bw)) local.bw <- NULL
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  cv.block <- floor(abs(cv.block))
  y <- as.matrix(y)
  x <- as.matrix(x)
  neq <- NCOL(y)
  obs <- NROW(y)
  if( !identical(obs, NROW(x)))
    stop("The number of equations in 'x' and 'y' are different \n")
  if(is.null(z)) 
  {
    lower <- 5 / obs
    upper <- ifelse(tkernel == "Gaussian", 10, 2.5)
  } 
  else 
  {
    z_range <- max(z) - min(z)
    z_sd <- stats::sd(z)
    lower <- max(z_range * 0.001, z_sd * 0.1) 
    upper <- z_range * ifelse(tkernel == "Gaussian", 10, 2.5)
  }
  weights <- NULL
  if (lmcv) 
  {
    z_grid <- if (is.null(z)) (1:obs)/obs else z
    center <- stats::median(z_grid)
    if (is.null(local.bw) && is.null(z)) 
    {
      local.bw <- 0.2 * (max(z_grid) - min(z_grid))
    }
    else if (is.null(local.bw) && !is.null(z))
    {
      z_iqr <- stats::IQR(z) / 1.349
      z_scale <- if (z_iqr > 0) min(stats::sd(z), z_iqr) else stats::sd(z)
      local.bw <- 0.9 * z_scale * (obs)^(-1/5)
    }
    u <- (z_grid - center)/local.bw
    if (tkernel == "Gaussian") 
    {
      weights <- stats::dnorm(u)
    } 
    else if (tkernel == "Epa") 
    {
      weights <- 0.75 * (1 - u^2) * (abs(u) <= 1)
    } 
    else if (tkernel == "Triweight") 
    {
      weights <- (35/32) * (1 - u^2)^3 * (abs(u) <= 1)
    }
    sum_w <- sum(weights)
    if (is.na(sum_w) || sum_w < .Machine$double.eps) 
    {
      warning("Sparse region, back to leave-one-out CV.")
      weights <- rep(1/obs, obs)
    } 
    else 
    {
      weights <- weights/sum_w
    }
  }
  find_grid <- function(h) 
  {
    val <- try(.tvOLS.cv(bw = h, x = x, y = y[, j], z = z, 
                         cv.block = cv.block, est = est, 
                         tkernel = tkernel, weights = weights, 
                         singular.ok = singular.ok), 
               silent = TRUE)
    if(inherits(val, "try-error") || is.na(val)) return(.Machine$double.xmax)
    return(val)
  }
  grid.starts <- exp(seq(log(lower), log(upper), length.out = 7))
  bw <- numeric(neq)
  for (j in 1:neq)
  {
    cv.values <- sapply(grid.starts, find_grid)
    best.start <- grid.starts[which.min(cv.values)]
    result <- try(stats::optim(par = best.start, fn = .tvOLS.cv, method = "Brent",
                               lower = lower, upper = upper, x = x, y = y[, j], z = z,
                               cv.block = cv.block, est = est, tkernel = tkernel,
                               weights = weights, singular.ok = singular.ok),
                  silent = TRUE)
    if (inherits(result, "list") && !is.na(result$value)) 
    {
      bw[j] <- result$par
    } 
    else 
    {
      bw[j] <- upper
      warning("Bandwidth optimization failed for equation ", j, ". Using classical OLS.")
    }
    if (isTRUE(all.equal(bw[j], lower))) 
      warning("Optimal bandwidth for equation ", j, " is at the lower bound. Results may be undersmoothed.")
  }
  return(bw)
}
#' @examples
#' data( Kmenta, package = "systemfit" )
#' 
#' ## x is a list of matrices containing the regressors, one matrix for each equation
#' x <- list()
#' x[[1]] <- Kmenta[, c("price", "income")]
#' x[[2]] <- Kmenta[, c("price", "farmPrice", "trend")]
#'
#' ## 'y' is a matrix with one column for each equation
#' y <- cbind(Kmenta$consump, Kmenta$consump)
#'
#' ## Select bandwidth by cross-validation
#' bw <- bw(x = x, y = y)
#'
#' ##One bandwidth per equation
#' print(bw)
#'
#' @rdname bw
#' @method bw list
#' @export
bw.list <- function(x, y, z = NULL, cv.block = 0, lmcv = FALSE,
                    local.bw = NULL, est = c("lc", "ll"), 
                    tkernel = c("Triweight", "Epa", "Gaussian"), 
                    singular.ok = TRUE, ...)
{
  if(!inherits(x, "list"))
    stop("'x' should be a list of matrices. \n")
  neq <- length(x)
  if(neq < 2)
    stop("'x' should have at least two elements.\n")
  if(sum(is.na(y)) > 0)
    stop("There are NA values in your data, please enter only complete cases. \n")
  if(is.null(y))
    stop("Parameter 'y' is missing.\n")
  obs <- NROW(x[[1]])
  if(!identical(neq, NCOL(y)) | !identical(obs, NROW(y)))
    stop("The number of equations in 'x' and 'y' are different \n")
  if (is.logical(local.bw)) local.bw <- NULL
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  cv.block <- floor(abs(cv.block))
  if(is.null(z)) 
  {
    lower <- 5 / obs
    upper <- ifelse(tkernel == "Gaussian", 10, 2.5)
  } 
  else 
  {
    z_range <- max(z) - min(z)
    z_sd <- stats::sd(z)
    lower <- max(z_range * 0.001, z_sd * 0.1) 
    upper <- z_range * ifelse(tkernel == "Gaussian", 10, 2.5)
  }
  weights <- NULL
  if (lmcv) 
  {
    z_grid <- if (is.null(z)) (1:obs)/obs else z
    center <- stats::median(z_grid)
    if (is.null(local.bw) && is.null(z)) 
    {
      local.bw <- 0.2 * (max(z_grid) - min(z_grid))
    }
    else if (is.null(local.bw) && !is.null(z))
    {
      z_iqr <- stats::IQR(z) / 1.349
      z_scale <- if (z_iqr > 0) min(stats::sd(z), z_iqr) else stats::sd(z)
      local.bw <- 0.9 * z_scale * (obs)^(-1/5)
    }
    u <- (z_grid - center)/local.bw
    if (tkernel == "Gaussian") 
    {
      weights <- stats::dnorm(u)
    } 
    else if (tkernel == "Epa") 
    {
      weights <- 0.75 * (1 - u^2) * (abs(u) <= 1)
    } 
    else if (tkernel == "Triweight") 
    {
      weights <- (35/32) * (1 - u^2)^3 * (abs(u) <= 1)
    }
    sum_w <- sum(weights)
    if (is.na(sum_w) || sum_w < .Machine$double.eps) 
    {
      warning("Sparse region, back to leave-one-out CV.")
      weights <- rep(1/obs, obs)
    } 
    else 
    {
      weights <- weights/sum_w
    }
  }
  
  find_grid <- function(h) 
  {
    val <- try(.tvOLS.cv(bw = h, x = x[[j]], y = y[, j], z = z, 
                         cv.block = cv.block, est = est, 
                         tkernel = tkernel, weights = weights, # Passed weights
                         singular.ok = singular.ok), 
               silent = TRUE)
    if(inherits(val, "try-error") || is.na(val)) return(.Machine$double.xmax)
    return(val)
  }
  grid.starts <- exp(seq(log(lower), log(upper), length.out = 7))
  bw <- numeric(neq)
  for (j in 1:neq)
  {
    cv.values <- sapply(grid.starts, find_grid)
    best.start <- grid.starts[which.min(cv.values)]
    result <- try(stats::optim(par = best.start, fn = .tvOLS.cv, 
                               method = "Brent", lower = lower, 
                               upper = upper, x = x[[j]], y = y[, j], z = z,
                               cv.block = cv.block, est = est, tkernel = tkernel, 
                               weights = weights, # Passed weights
                               singular.ok = singular.ok), silent = TRUE)
    if (inherits(result, "list") && !is.na(result$value)) 
    {
      bw[j] <- result$par
    } 
    else 
    {
      bw[j] <- upper
      warning("Bandwidth optimization failed for equation ", j, ". Using Global OLS.")
    }
    if (isTRUE(all.equal(bw[j], lower))) 
      warning("Optimal bandwidth for equation ", j, " is at the lower bound. Results may be undersmoothed.")
  }
  return(abs(bw))
}

#' @rdname bw
#' @method bw tvlm
#' @export
bw.tvlm <- function(x, ...)
{
  y <- x$y
  z <- x$z
  cv.block <- floor(abs(x$cv.block))
  lmcv <- x$lmcv
  local.bw <- x$local.bw
  est <- x$est
  tkernel <- x$tkernel
  singular.ok <- x$singular.ok
  x <- x$x
  return(bw (x, y, z, cv.block, lmcv, local.bw, est, tkernel, singular.ok))
}

#' @rdname bw
#' @method bw tvar
#' @export
bw.tvar <- bw.tvlm

#' @rdname bw
#' @method bw tvvar
#' @export
bw.tvvar <- bw.tvlm

#' @rdname bw
#' @method bw tvsure
#' @export
bw.tvsure <- bw.tvlm

#' @rdname bw
#' @method bw tvplm
#' @export
bw.tvplm <- function(x, ...)
{
  if(!inherits(x, "tvplm"))
    stop("'x' should be a 'tvplm' object. \n")
  y <- x$y
  z <- x$z
  cv.block <- floor(abs(x$cv.block))
  obs <- x$obs
  neq <- x$neq
  method <- x$method
  est <- x$est
  tkernel <- x$tkernel
  x <- x$x
  if(is.null(z))
  {
    upper <- 2.5
    lower <- 5/obs
    top <- 1
  }
  else
  {
    top <- (max(z) - min(z))*5
    upper <- top
    lower <- top * 0.001
  }
  iter <- 0
  value <- .Machine$double.xmax
  bw <- top
  while(value == .Machine$double.xmax)
  {
    if(iter == 10)
    {
      warning("Maximum number of iterations reached in bandwidth calculation: either the function
              is constant, or no convergence of bandwidth. \n")
      bw <- top
      break()
    }
    if (method != "within")
      result <- try(stats::optim(stats::runif(1, lower, top), .tvRE.cv, method = "Brent",
                                 lower = lower, upper = upper, x = x, y = y, z = z, 
                                 neq = neq, obs = obs, cv.block = cv.block, est = est, 
                                 tkernel = tkernel),
                    silent = TRUE)
    else 
      result <- try(stats::optim(stats::runif(1, lower, top), .tvFE.cv, method = "Brent",
                                 lower = lower, upper = top, x = x, y = y, z = z, 
                                 neq = neq, obs = obs, cv.block = cv.block,  est = est, 
                                 tkernel = tkernel),
                    silent = TRUE)
    if (inherits(result, "list") && !is.na(result$value))
    {
      value <- result$value 
      bw <- result$par
    } 
    else 
    {
      value <- .Machine$double.xmax
    }
    iter <- iter + 1
  }
  return(abs(bw))  
}


#' @rdname bw
#' @method bw pdata.frame
#' @export
bw.pdata.frame <- function(x, z = NULL, method, cv.block = 0, 
                           est = c("lc", "ll"), 
                           tkernel = c("Triweight", "Epa", "Gaussian"), ...)
{  
  dimen <- plm::pdim(x)
  neq <- dimen$nT$n
  obs <- dimen$nT$T
  y <- stats::model.extract(x, "response")
  if (stats::is.empty.model(x))
    stop ("No regressors in the model. \n")
  terms <- attr(x, "terms")
  x_mat <- stats::model.matrix(terms, x)
  
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  if(is.null(z)) 
  {
    lower <- 5 / obs
    upper <- ifelse(tkernel == "Gaussian", 10, 2.5)
  } 
  else 
  {
    z_range <- max(z) - min(z)
    z_sd <- stats::sd(z)
    lower <- max(z_range * 0.001, z_sd * 0.1) 
    upper <- z_range * ifelse(tkernel == "Gaussian", 10, 2.5)
  }
  cv.func <- if (method != "within") .tvRE.cv else .tvFE.cv
  
  find_grid <- function(h) 
  {
    val <- try(cv.func(bw = h, x = x_mat, y = y, z = z, 
                       neq = neq, obs = obs, cv.block = cv.block, 
                       est = est, tkernel = tkernel), 
               silent = TRUE)
    if(inherits(val, "try-error") || is.na(val)) return(.Machine$double.xmax)
    return(val)
  }
  grid.starts <- exp(seq(log(lower), log(upper), length.out = 7))
  cv.values <- sapply(grid.starts, find_grid)
  best.start <- grid.starts[which.min(cv.values)]  
  result <- try(stats::optim(par = best.start, fn = cv.func, method = "Brent",
                                 lower = lower, upper = upper, x = x_mat, y = y, z = z, 
                                 neq = neq, obs = obs, cv.block = cv.block, est = est, 
                                 tkernel = tkernel),
                    silent = TRUE)
  if (inherits(result, "list") && !is.na(result$value)) 
  {
    bw <- result$par
  } 
  else 
  {
    bw <- upper
    warning("Bandwidth optimization failed for pdata.frame model. Using Global OLS (upper bound).")
  }
  if (isTRUE(all.equal(bw, lower))) 
    warning("Optimal bandwidth is at the lower bound. Results may be undersmoothed.")
  return(abs(bw))   
}


#' Covariance Bandwidth Calculation by Cross-Validation
#' 
#' \code{bwCov} calculates a single bandwidth to estimate the time-varying variance-
#' covariance matrix.
#' 
#' @param x A matrix or a data frame.
#' @inheritParams bw
#' @return A scalar.
#' @examples
#'
#' data(CEES)
#' ## Using a shorter set for a quick example. Variable "Date" is removed.
#' mydata <- tail (CEES[, -1], 50)
#' bw.cov <- bwCov(mydata)
#' Sigma.hat <- tvCov(mydata, bw = bw.cov)
#'
#' @rdname bwCov
#' @export
bwCov <- function(x, z = NULL, cv.block = 0, est = c("lc", "ll"), 
                  tkernel = c("Triweight", "Epa", "Gaussian"))
{
  if(!inherits(x, c("matrix", "data.frame")))
    stop("'x' should be a matrix or a data.frame.\n")
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  x <- as.matrix(x)
  cv.block <- abs(cv.block)
  obs <- NROW(x)
  neq <- NCOL(x)
  if(is.null(z)) 
  {
    lower <- 5 / obs
    upper <- ifelse(tkernel == "Gaussian", 10, 2.5)
  } 
  else 
  {
    z_range <- max(z) - min(z)
    z_sd <- stats::sd(z)
    lower <- max(z_range * 0.001, z_sd * 0.1) 
    upper <- z_range * ifelse(tkernel == "Gaussian", 10, 2.5)
  }
  find_grid <- function(h) 
  {
    val <- try(.tvCov.cv(bw = h, x = x, z = z, cv.block = cv.block, 
                         est = est, tkernel = tkernel), 
               silent = TRUE)
    if(inherits(val, "try-error") || is.na(val)) return(.Machine$double.xmax)
    return(val)
  }
  grid.starts <- exp(seq(log(lower), log(upper), length.out = 7))
  cv.values <- sapply(grid.starts, find_grid)
  best.start <- grid.starts[which.min(cv.values)]
  result <- try(stats::optim(par = best.start, fn = .tvCov.cv, method = "Brent",
                              lower = lower, upper = upper, x = x, z = z, cv.block = cv.block,
                              est = est, tkernel = tkernel),
                  silent = TRUE)
  if (inherits(result, "list") && !is.na(result$value)) 
  {
    bw <- result$par
  } 
  else 
  {
    bw <- upper
    warning("Bandwidth optimization failed for Covariance. Using Global Covariance (upper bound).")
  }
  
  if (isTRUE(all.equal(bw, lower))) 
    warning("Optimal bandwidth for Covariance is at the lower bound. Results may be undersmoothed.")
  return(abs(bw))
}