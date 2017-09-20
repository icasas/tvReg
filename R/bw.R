##' Bandwidth Selection by Cross-Validation
##'
##' Calculate bandwidth(s) by cross-validation for functions tvSURE, tvVAR and tvLM.
##'
##' @rdname bw
##' @param x an object used to select a method.
##' @param ... Other parameters passed to specific methods.
##' @return \code{bw} returns a vector or a scalar with the bandwith to estimate the mean or the covariance
##' residuals, fitted values.
##' @export
bw <- function(x, ...)  UseMethod("bw", x)


##' @rdname bw
##' @param y A matrix or vector with the dependent variable(s).
##' @param z A vector with the variable over which coefficients are smooth over.
##' @param est The nonparametric estimation method, one of "lc" (default) for linear constant
##' or "ll" for local linear.
##' @param tkernel The type of kernel used in the coefficients estimation method,
##' one of Epanesnikov ("Epa") or "Gaussian".
##' @param singular.ok	Logical. If FALSE, a singular model is an error.
##'
##' @return A scalar or a vector of scalars.
##'
##' @examples
##' tau <- seq(1:1000)/1000
##' beta <- data.frame(beta1 = sin(2*pi*tau), beta2 =  2*tau)
##' X <- data.frame(X1 = rnorm(1000), X2 =  rchisq(1000, df = 4))
##' error <- rt(1000, df = 10)
##' y <- apply(X*beta, 1, sum) + error
##' bw <- bw(X, y, est = "ll")
##' model.tv <-  tvOLS(x = X, y = y, bw = bw)
##'
##' @method bw default
##' @export
##'
bw.default <- function(x, y, z = NULL, tkernel = "Epa", est = "lc",
                       singular.ok = TRUE, ...)
{
  if(!is.matrix(x) & !is.data.frame(x) & !is.numeric(x) & !is.vector(x))
    stop("'x' should be a matrix, a vector or a data frame. \n")
  if(is.null(y))
    stop("Parameter 'y' missing. \n")
  if(sum(is.na(y))>0 | sum(is.na(x))>0)
    stop("There are NA values in your data, please enter only complete cases. \n")
  y <- as.matrix(y)
  x <- as.matrix(x)
  neq <- ncol(y)
  obs <- nrow(y)
  if(is.null(z))
  {
    upper <- 20
    lower <- 5/obs
    top <- 1
  }
  else
  {
    top <- max(z)- min(z)
    upper <- top * 5
    dist <- diff(sort(z))
    lower <- min(dist)
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
        break()
      }
      result <- try(stats::optim(stats::runif(1, lower, top), .tvOLS.cv, method = "Brent",
                                 lower = lower, upper = upper, x = x, y = y[, j], z = z,
                                 est = est, tkernel = tkernel, singular.ok = singular.ok),
                    silent = TRUE)
      if (class(result)  !=  "list")
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
  if(iter == 10)
    warning("Maximum number of iterations reached in bandwidth calculation: either the function
            is constant, or no convergence of bandwidth. \n")
  return(bw)
}


##' @examples
##' data( Kmenta, package = "systemfit" )
##' ## x is a list of matrices containing the regressors, one matrix for each equation
##' x <- list()
##' x[[1]] <- Kmenta[, c("price", "income")]
##' x[[2]] <- Kmenta[, c("price", "farmPrice", "trend")]
##'
##' ## y is a matrix with one column for each equation
##' y <- cbind(Kmenta$consump, Kmenta$consump)
##'
##' ## Calculate bandwidth
##' bw <- bw(x = x, y = y, tkernel = "Gaussian")
##'
##' ##One bandwidth per equation
##' print(bw)
##'
##' ##Use these bandwidths to estimate the time-varying coefficients
##' tvgls <- tvGLS(x = x, y = y, bw = bw)
##'
##' @rdname bw
##' @method bw list
##' @export
bw.list <- function(x, y, z = NULL, est = "lc", tkernel = "Epa",
                    singular.ok = singular.ok, ...)
{
  if(!is.list(x))
    stop("'x' should be a list of matrices. \n")
  neq <- length(x)
  if(neq< 2)
    stop("'x' should have at least two elements.\n")
  if (sum(is.na(y))>0)
    stop("There are NA values in your data, please enter only complete cases. \n")
  if(is.null(y))
    stop("Parameter 'y' is missing.\n")
  obs <- nrow(x[[1]])
  if(!identical(neq, length(x)[1]))
    stop("The number of equations in 'x' and 'y' are different \n")
  if (!identical(obs, nrow(x[[1]])))
    stop("The length of 'x' and 'y' are different \n")
  if (est  !=  "lc" & est  !=  "ll")
    stop("The supported estimation methods in parameter 'est' are: 'lc' and 'll' \n")
  if (tkernel  !=  "Epa" & tkernel  !=  "Gaussian")
    stop("The supported kernels are 'Epa' and 'Gaussian'\n")
  bw <- numeric(neq)
  if(is.null(z))
  {
    upper <- 2
    lower <- 5/obs
    top <- 0.5
  }
  else
  {
    top <- max(z)- min(z)
    upper <- top * 5
    dist <- diff(sort(z))
    lower <- min(dist)
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
        break()
      }
      result <- try(stats::optim(stats::runif(1, lower, top), .tvOLS.cv, method = "Brent",
                                 lower = lower, upper = upper, x = x[[j]], y = y[, j],
                                 est = est, tkernel = tkernel, singular.ok = singular.ok),
                    silent = TRUE)
      if (class(result)  !=  "list")
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
  if(iter == 10)
    warning("Maximum number of iterations reached in bandwidth calculation: either the function
            is constant, or no convergence of bandwidth. \n")
  return(bw)
}


##' Covariance Bandwidth Calculation by Cross-Validation
##' \emph{bwCov} calculates a single bandwidth to estimate the time-varying variance-
##' covariance matrix.
##'
##' @param x A matrix or a data frame.
##' @inheritParams bw
##' @return A scalar.
##' @examples
##'
##' data(CEES)
##' ## Using a shorter set for a quick example
##' mydata <- tail (CEES, 100)
##' bw.cov <- bwCov(mydata)
##' Sigma.hat <-  tvCov(mydata, bw = bw.cov)
##'
##' @rdname bwCov
##' @export
bwCov <- function(x, tkernel = "Epa",  est = "lc")
{
  if(!is.matrix(x) & !is.data.frame(x))
    stop("'x' should be a matrix or a data.frame.\n")
  x <- as.matrix(x)
  obs <- nrow(x)
  neq <- ncol(x)
  value  <- .Machine$double.xmax
  iter <- 0
  while(value == .Machine$double.xmax)
  {
    if(iter == 10)
    {
      value <- 0
      bw <- 100
      break()
    }
    result <- try(stats::optim(stats::runif(1, 5/obs, 1), .tvCov.cv, method = "Brent",
                               lower = 5/obs, upper = 20, x = x, est = est,
                               tkernel = tkernel),
                  silent = TRUE)
    if (class(result)  !=  "list")
      value <- .Machine$double.xmax
    else
    {
      value <- result$value
      bw <- result$par
    }
    iter <- iter + 1
  }
  if (iter == 10)
    warning("Maximum number of iterations reached in bandwidth calculation: either the function
            is constant, or no convergence of bandwidth. \n")
  return(bw)
}

