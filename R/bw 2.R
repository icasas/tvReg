#' Bandwidth Selection by Cross-Validation
#'
#' Calculate bandwidth(s) by cross-validation for functions tvSURE, tvVAR and tvLM.
#'
#' @rdname bw
#' @param x An object used to select a method.
#' @param ... Other parameters passed to specific methods.
#' @return \code{bw} returns a vector or a scalar with the bandwith to estimate the mean or the covariance
#' residuals, fitted values.
#' @export
bw <- function(x, ...)  UseMethod("bw", x)


#' @rdname bw
#' @param y A matrix or vector with the dependent variable(s).
#' @param z A vector with the variable over which coefficients are smooth over.
#' @param cv.block A positive scalar with the size of the block in leave-one block-out cross-validation.
#' By default 'cv.block=0' meaning leave-one-out cross-validation.
#' @param est The nonparametric estimation method, one of "lc" (default) for linear constant
#' or "ll" for local linear.
#' @param tkernel A character, either "Triweight" (default), "Epa" or "Gaussian" kernel function.
#' @param singular.ok	Logical. If FALSE, a singular model is an error.
#'
#' @return A scalar or a vector of scalars.
#'
#' @examples
#' ##Generate data
#' tau <- seq(1:200)/200
#' beta <- data.frame(beta1 = sin(2*pi*tau), beta2 =  2*tau)
#' X <- data.frame(X1 = rnorm(200), X2 =  rchisq(200, df = 4))
#' error <- rt(200, df = 10)
#' y <- apply(X*beta, 1, sum) + error
#' 
#' ##Select bandwidth by cross-validation
#' bw <- bw(X, y, est = "ll", tkernel = "Gaussian")
#'
#' @method bw default
#' @export
#'
bw.default <- function(x, y, z = NULL, cv.block = 0, est = c("lc", "ll"), tkernel = c("Triweight", "Epa", "Gaussian"),
                       singular.ok = TRUE, ...)
{ 
  if(!inherits(x, c("matrix", "data.frame", "vector", "numeric", "integer")))
    stop("'x' should be a matrix, a vector or a data frame. \n")
  if(is.null(y))
    stop("Parameter 'y' missing. \n")
  if(sum(is.na(y))>0 | sum(is.na(x))>0)
    stop("There are NA values in your data, please enter only complete cases. \n")
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
    upper <- 20
    lower <- 5/obs
    top <- 1
  }
  else
  {
    top <- (max(z) - min (z))* 5
    upper <- top 
    lower <- top * 0.001
  }
  bw <- numeric(neq)
  for (j in 1:neq)
  {
    iter <- 0
    value <- .Machine$double.xmax
    while(value == .Machine$double.xmax)
    {
      if(iter == 10)
      {
        value <- 0
        bw[j] <- top
        warning("Maximum number of iterations reached in bandwidth calculation: either the function
            is constant, no convergence of bandwidth, or cv.block is too big. \n")
        break()
      }
      result <- try(stats::optim(stats::runif(1, lower, top), .tvOLS.cv, method = "Brent",
                                 lower = lower, upper = upper, x = x, y = y[, j], z = z,
                                 cv.block = cv.block, est = est, tkernel = tkernel, 
                                 singular.ok = singular.ok),
                    silent = TRUE)
      if (!inherits(result, "list"))
        value <- .Machine$double.xmax
      else
      {
        if(is.na(result$value))
          value <- .Machine$double.xmax
        else
        {
          value <- result$value
          bw[j] <- result$par
        }
      }
      iter <- iter + 1
    }
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
bw.list <- function(x, y, z = NULL, cv.block = 0, est = c("lc", "ll"), tkernel = c("Triweight", "Epa", "Gaussian"),
                    singular.ok = TRUE, ...)
{
  if(!inherits(x, "list"))
    stop("'x' should be a list of matrices. \n")
  neq <- length(x)
  if(neq < 2)
    stop("'x' should have at least two elements.\n")
  if (sum(is.na(y)) >0)
    stop("There are NA values in your data, please enter only complete cases. \n")
  if(is.null(y))
    stop("Parameter 'y' is missing.\n")
  obs <- NROW(x[[1]])
  if(!identical(neq, NCOL(y)) | !identical(obs, NROW(y)))
    stop("The number of equations in 'x' and 'y' are different \n")
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  bw <- numeric(neq)
  cv.block <- floor(abs(cv.block))
  if(is.null(z))
  {
    upper <- 20
    lower <- 5/obs
    top <- 0.5
  }
  else
  {
    top <- (max(z) - min(z))*5
    upper <- top
    lower <- top * 0.001
  }
  for (j in 1:neq)
  {
    iter <- 0
    value <- .Machine$double.xmax
    while(value == .Machine$double.xmax)
    {
      if(iter == 10)
      {
        value <- 0
        bw[j] <- 100
        warning("Maximum number of iterations reached in bandwidth calculation: either the function
            is constant, or no convergence of bandwidth. \n")
        break()
      }
      result <- try(stats::optim(stats::runif(1, lower, top), .tvOLS.cv, 
                                 method = "Brent", lower = lower, 
                                 upper = upper, x = x[[j]], y = y[, j], z = z,
                                 cv.block = cv.block, est = est, tkernel = tkernel, 
                                 singular.ok = singular.ok), silent = TRUE)
      if (!inherits(result, "list"))
        value <- .Machine$double.xmax
      else
      {
        if(is.na(result$value))
          value <- .Machine$double.xmax
        else
        {
          value <- result$value
          bw[j] <- result$par
        }
      }
      iter <- iter + 1
    }
  }
  return(bw)
}
#' @rdname bw
#' @method bw tvlm
#' @export

bw.tvlm <- function(x, ...)
{
  y <- x$y
  z <- x$z
  est <- x$est
  tkernel <- x$tkernel
  singular.ok <- x$singular.ok
  cv.block <- floor(abs(x$cv.block))
  x <- x$x
  return(bw (x, y, z, cv.block, est, tkernel, singular.ok))
}

#' @rdname bw
#' @method bw tvar
#'
#' @export
bw.tvar <- bw.tvlm

#' @rdname bw
#' @method bw tvvar
#'
#' @export
bw.tvvar <- bw.tvlm

#' @rdname bw
#' @method bw tvsure
#'
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
    upper <- 20
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
  while(value == .Machine$double.xmax)
  {
    if(iter == 10)
    {
      value <- 0
      bw <- 100
      warning("Maximum number of iterations reached in bandwidth calculation: either the function
              is constant, or no convergence of bandwidth. \n")
      break()
    }
    if (method != "within")
      result <- try(stats::optim(stats::runif(1, lower, top), .tvRE.cv, method = "Brent",
                                 lower = lower, upper = upper, x = x, y = y, z = z, 
                                 neq = neq, obs = obs, cv.block = cv.block, est = est, 
                                 tkernel = tkernel),
                    silent = FALSE)
    else 
      result <- try(stats::optim(stats::runif(1, lower, top), .tvFE.cv, method = "Brent",
                                 lower = lower, upper = top, x = x, y = y, z = z, 
                                 neq = neq, obs = obs, cv.block = cv.block,  est = est, 
                                 tkernel = tkernel),
                    silent = FALSE)
    if (!inherits(result, "list"))
      value <- .Machine$double.xmax
    else
    {
      if(is.na(result$value))
        value <- .Machine$double.xmax
      else
      {
        value <- result$value 
        bw <- result$par
      }
    }
    iter <- iter + 1
  }
  return(abs(bw))  
}

#' 
#' Panel Model Bandwidth Calculation by Cross-Validation
#' \emph{bwPanel} calculates a single bandwidth to estimate the time-varying 
#' coefficients of a panel data modelc
#' @param method A character with the choice of panel model/estimation method:
#' If method = \code{tvPOLS} (default) then the data is pooled estimated with time-varying OLS. 
#' No individual or time effects are estimated
#' If method = \code{tvFE} then individual effects which might be correlated with 
#' the regressors are estimated.
#' If method = \code{tvRE} then individual effects are considered random and independent
#' of the regressors.
#' @return A scalar.
#' @method bw pdata.frame
#' @rdname bw
#' @export
bw.pdata.frame<-function(x, z = NULL, method, cv.block = 0, 
                  est = c("lc", "ll"), tkernel = c("Triweight", "Epa", "Gaussian"), ...)
{  
  dimen <- plm::pdim(x)
  neq <- dimen$nT$n
  obs <- dimen$nT$T
  y <- stats::model.extract(x, "response")
  if (stats::is.empty.model(x))
    stop ("No regressors in the model. \n")
  else
  {
    terms <- attr(x, "terms")
    x <- stats::model.matrix(terms, x)
    var.names <- colnames(x)
  }
  nvar <-  NCOL (x)
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  if(is.null(z))
  {
    upper <- 20
    lower <- 5/obs
    top <- 1
  }
  else
  {
    top <- (max(z) - min(z))*5
    upper <- top
    lower <- top * 0.001
  }
  value <- .Machine$double.xmax
  iter <- 0
  while(value == .Machine$double.xmax)
  {
    if(iter == 10)
    {
      value <- 0
      bw <- top
      warning("Maximum number of iterations reached in bandwidth calculation: either the function
            is constant, or no convergence of bandwidth. \n")
      break()
    }

    if (method != "within")
      result <- try(stats::optim(stats::runif(1, lower, top), .tvRE.cv, method = "Brent",
                                 lower = lower, upper = upper, x = x, y = y, z = z, 
                                 neq = neq, obs = obs, cv.block = cv.block, est = est, 
                                 tkernel = tkernel),
                    silent = FALSE)
    else 
      result <- try(stats::optim(stats::runif(1, lower, top), .tvFE.cv, method = "Brent",
                                 lower = lower, upper = top, x = x, y = y, z = z, 
                                 neq = neq, obs = obs, cv.block = cv.block,  est = est, 
                                 tkernel = tkernel),
                    silent = FALSE)
    if (!inherits(result, "list"))
        value <- .Machine$double.xmax
    else
    {
      if(is.na(result$value))
        value <- .Machine$double.xmax
      else
      {
        value <- result$value 
        bw <- result$par
      }
    }
    iter <- iter + 1
  }
  return(abs(result$par))  
}



#' 
#' Covariance Bandwidth Calculation by Cross-Validation
#' \emph{bwCov} calculates a single bandwidth to estimate the time-varying variance-
#' covariance matrix.
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
bwCov <- function(x, cv.block = 0, est = c("lc", "ll"), tkernel = c("Triweight", "Epa", "Gaussian"))
{
  if(!inherits(x, c("matrix", "data.frame")))
    stop("'x' should be a matrix or a data.frame.\n")
  tkernel <- match.arg(tkernel)
  est <- match.arg(est)
  x <- as.matrix(x)
  cv.block <- abs(cv.block)
  obs <- NROW(x)
  neq <- NCOL(x)
  value  <- .Machine$double.xmax
  iter <- 0
  while(value == .Machine$double.xmax)
  {
    if(iter == 10)
    {
      value <- 0
      bw <- 20
      warning("Maximum number of iterations reached in bandwidth calculation: either the function
            is constant, or no convergence of bandwidth. \n")
      break()
    }
    result <- try(stats::optim(stats::runif(1, 5/obs, 1), .tvCov.cv, method = "Brent",
                               lower = 5/obs, upper = 20, x = x, cv.block = cv.block,
                               est = est, tkernel = tkernel),
                  silent = TRUE)
    if (!inherits(result, "list"))
      value <- .Machine$double.xmax
    else
    {
      value <- result$value
      bw <- result$par
    }
    iter <- iter + 1
  }
  return(bw)
}
