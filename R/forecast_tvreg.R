
#' Forecast Methods for Objects in tvReg.
#'
#' \code{forecast} calculates the forecast for objects with class attribute \code{tvlm}, \code{tvar}, 
#' \code{tvvar}, \code{tvirf}, \code{tvsure} and \code{tvplm}. If the 
#' smoothing variable (z) in the model is non-NULL and it is a random 
#' variable then use function \code{predict} with parameter \code{newz}.
#' @param object An object used to select a method.
#' @param ... Other parameters passed to specific methods.
#' @return An object of class matrix or vector with the same dimensions than  the dependent 
#' variable of \code{object}.
#' @seealso \code{\link{predict}}. 
#' @rdname forecast-tvReg
#' @import methods
#' @export
forecast <- function(object, ...) UseMethod("forecast", object)

#' @rdname forecast-tvReg
#' @method forecast tvlm
#' @param newdata A data.frame, matrix or vector with new values of all independent variables. No need to 
#' input the intercept. Note that this does not refer to the variables in 'exogen' which might be part of 'tvar' 
#' and 'tvvar' objects. Those must be included in 'newexogen'.
#' @param n.ahead A scalar with the forecast horizon, value 1 by default.
#' @param winsize A scalar. If 0 then an 'increase window' forecasting is performed.
#'  Otherwise a 'rolling window' forecasting is performed with window size given by 
#'  'winsize'.
#' @examples 
#' data("RV")
#' RV2 <- head(RV, 2001)
#' TVHAR <- tvLM (RV ~ RV_lag + RV_week + RV_month, data = RV2, bw = 20)
#' newdata <- cbind(RV$RV_lag[2002:2004], RV$RV_week[2002:2004],
#'               RV$RV_month[2002:2004])
#' forecast(TVHAR, newdata, n.ahead = 3)
#' 
#' @export
forecast.tvlm<-function (object, newdata, n.ahead = 1, winsize = 0, ...) 
{
  if(!inherits(object, c("tvlm")))
    stop("\nParameter 'object' should be entered and it should have class 'tvlm' or 'tvar'.\n")
  if(!is.null(object$z))
    stop("\nYour model coefficients are functions of a random variable 'z', use function 
         'predict' with parameter 'newz'.\n")
  if(!inherits(newdata, c("data.frame", "matrix", "numeric", "vector")))
    stop("\nParameter 'newdata' should be entered and it should be a numeric vector if there is only
         one row or a 'matrix' or a 'data.frame' for more than one row.\n")
  if(n.ahead == 1) 
    newdata <- matrix(newdata, ncol = length(newdata))
  if(NROW(newdata) != n.ahead)
    stop("\nDimensions of 'newdata' are not compatible with 'n.ahead'.\n")
  n.col <- NCOL(newdata)
  is.intercept <- ("(Intercept)" %in% colnames(object$x))
  if(is.intercept & n.col ==(NCOL(object$x) - 1))
    newdata <- cbind(rep(1, n.ahead), newdata)
  obs <- object$obs
  is.rw <- !(winsize == 0)
  winsize <- abs (winsize)
  if(winsize > obs)
    winsize <- obs - 1
  prediction <- numeric(n.ahead)
  totobs <- obs + n.ahead
  grid <- (1:totobs)/totobs
  X <- object$x
  Y <- object$y
  for (t in 1:n.ahead)
  {
    i <- ifelse (is.rw, obs - winsize, 1)
    object$x <- X[i:(obs + t -1),]
    object$y <- Y[i:(obs + t -1)]
    object$z <- grid[i:(obs + t -1)]
    object$ez <- grid[obs + t]
    theta <- try(tvOLS(object)$coefficients)
    if (inherits(theta, "try-error"))
    {
      object$bw <- bw (object)
      theta<- tvOLS(object)$coefficients
    }
    prediction[t] <- theta%*%newdata[t,]
    X <- rbind(X, newdata[t,])
    Y <- c (Y, prediction[t])
  }
  return(prediction)
}

#' @rdname forecast-tvReg
#' @method forecast tvar
#' @param newz A vector with the new values of the smoothing variable.
#' @param newexogen A matrix or vector with the new values of the exogenous variables.
#' Only for predictions of *tvar* and *tvvar* objects.
#' @examples 
#' exogen = RV2[, c("RV_week", "RV_month")]
#' TVHAR2 <- tvAR(RV2$RV_lag, p = 1, exogen = exogen, bw = 20)
#' newexogen <- newdata[, -1]
#' forecast(TVHAR2, n.ahead = 3, newexogen = newexogen)
#' 
#' @export
forecast.tvar <- function(object, n.ahead = 1, newz = NULL, newexogen = NULL, winsize = 0, ...) 
{
  if(!inherits(object, c("tvar")))
    stop("\nParameter 'object' should be entered and it should have class 'tvar'.\n")
  if(!inherits(newz, c( "numeric", "vector")) & !is.null(object$z))
    stop("\nArgument 'newz' should be entered and it should be a numeric vector.\n")
  is.exogen <- !is.null(newexogen)
  if(is.exogen)
  {
    if (!identical(NCOL(newexogen), NCOL(object$exogen))) 
      stop("\nWrong dimension in 'newexogen'.\n")
    newexogen <- as.matrix(newexogen)
    if(NROW(newexogen) != n.ahead)
      stop("\nDimensions of 'newxexogen' and 'n.ahead' are not compatible.\n")
    newexogen <- as.matrix(newexogen)
  }
  is.intercept <- (object$type == "const")
  prediction <- numeric(n.ahead)
  obs <- object$obs
  is.rw <- !( winsize == 0)
  winsize <- abs (winsize)
  if(winsize > obs)
    winsize <- obs - 1
  prediction <- numeric(n.ahead)
  totobs <- obs + n.ahead
  if(is.null(newz))
    grid <- (1:totobs)/totobs
  else
    grid <- c(object$z, newz)
  X <- object$x
  Y <- object$y
  p <- object$p
  mask <- object$mask
  for (t in 1:n.ahead)
  {
    i <- ifelse (is.rw, obs - winsize, 1)
    object$x <- X[i:(obs + t -1),]
    object$y <- Y[i:(obs + t -1)]
    object$z <- grid[i:(obs + t -1)]
    object$ez <- grid[obs + t]
    theta <- try(tvOLS(object)$coefficients)
    if (inherits(theta, "try-error"))
    {
      object$bw <- bw (object)
      theta<- tvOLS(object)$coefficients
    }
    newx <- tail(object$y, p)
    if(is.intercept) newx <- c(1L, newx)
    if(is.exogen) newx <- c(newx, newexogen[t, ])
    newx <- newx[mask]
    prediction[t] <- theta%*%newx
    X <- rbind(X, newx)
    Y <- c (Y, prediction[t])
  }
  return(prediction)
}



#' @rdname forecast-tvReg
#' @method forecast tvvar
#' @examples 
#' data(usmacro, package = "bvarsv")
#' model.tvVAR <- tvVAR(usmacro, p = 6, type = "const", bw = c(1.8, 20, 20))
#' forecast(model.tvVAR, n.ahead = 10)
#' 
#' @export
forecast.tvvar<-function (object, n.ahead = 1, newz = NULL, newexogen = NULL, winsize = 0, ...) 
{
  if(!inherits(object, c("tvvar")))
    stop("\nParameter 'object' should be entered and it should have class 'tvvar'.\n")
  if(!inherits(newz, c( "numeric", "vector")) & !is.null(object$z))
    stop("\nArgument 'newz' should be entered and it should be a numeric vector.\n")
  if (!identical(NCOL(newexogen), NCOL(object$exogen))) 
    stop("\nWrong dimension in 'newexogen'.\n")
  is.intercept <- (object$type == "const")
  is.exogen <- FALSE
  if(!is.null(newexogen))
  {
    is.exogen <- TRUE
    if (!identical(NCOL(newexogen), NCOL(object$exogen))) 
      stop("\nWrong dimension in 'newexogen'.\n")
    if(!inherits(newexogen, c("numeric", "vector", "matrix", "data.frame")))
      stop("\nWrong class in 'newexogen'.\n")
    if(is.null(dim(newexogen)))
      newexogen <- matrix (newexogen, nrow = length(newexogen), ncol = 1)
    if(NROW(newexogen) != n.ahead)
      stop("\nDimensions of 'newxexogen' and 'n.ahead' are not compatible.\n")
  }
  neq <- object$neq
  obs <- object$obs
  prediction <- matrix(NA, nrow = n.ahead, ncol = neq)
  colnames(prediction) <- colnames(object$y)
  is.rw <- !(winsize == 0)
  winsize <- abs (winsize)
  if(winsize > obs)
    winsize <- obs - 1
  i <- ifelse (is.rw, obs - winsize, 1)
  totobs <- obs + n.ahead
  if(is.null(newz))
    grid <- (1:totobs)/totobs
  else
    grid <- c(object$z, newz)
  X <- object$x
  Y <- object$y
  p <- object$p
  nlags <- p*neq
  for (t in 1:n.ahead)
  {
    object$x <- X[i:(obs + t -1), ]
    object$y <- Y[i:(obs + t -1), ]
    object$z <- grid[i:(obs + t -1)]
    object$ez <- grid[obs + t]
    theta <- try(tvOLS(object)$coefficients, silent = TRUE)
    if (inherits(theta, "try-error"))
    {
      object$bw <- bw (object)
      theta<- tvOLS(object)$coefficients
    }
    rhs <- c(tail(object$y, 1))
    rhs <- c(rhs, tail(object$x, 1)[1:(nlags - neq)])
    if(is.intercept) 
      rhs <- c(rhs, "(Intercept)" = 1L)
    if(is.exogen)
      rhs <- c(rhs, newexogen[t, ])
    for (i in 1:neq)
      prediction[t, i] <- theta[[i]]%*%rhs
    X <- rbind(X, rhs)
    Y <- rbind(Y, prediction[t, ])
  }
  return(prediction)
}

#' @rdname forecast-tvReg
#' @method forecast tvsure
#' @param newdata A matrix or data.frame with the values of the regressors to use
#' for forecasting.
#' @examples 
#' data("Kmenta", package = "systemfit")
#' eqDemand <- consump ~ price + income
#' eqSupply <- consump ~ price + farmPrice 
#' system <- list(demand = eqDemand, supply = eqSupply)
#' tvOLS.fit <- tvSURE(system, data = Kmenta, est = "ll", bw = c(1.5, 1.5))
#' newdata <- data.frame(price = c(90, 100, 103), farmPrice = c(70, 95, 103), 
#' income = c(82, 94, 115))
#' forecast(tvOLS.fit, newdata = newdata, n.ahead = 3)
#' @export
forecast.tvsure<-function (object, newdata, n.ahead = 1, winsize = 0, ...) 
{
  if(!inherits(object, c("tvsure")))
    stop("\nParameter 'object' should be entered and it should have class 'tvsure'.\n")
  if(!is.null(object$z))
    stop("\nYour model coefficients are functions of a random variable 'z', use function 
         'predict' with parameter 'newz'.\n")
  if(!inherits(newdata, c("matrix", "data.frame")))
    stop("\nParameter 'newdata' should be a 'matrix' or a 'data.frame'.\n")
  if(is.null(colnames(newdata)))
    stop("\nParameter 'newdata' needs variable names.\n")
  if(NROW(newdata) != n.ahead)
    stop("\nDimensions of 'newdata' are not compatible with 'n.ahead'.\n")
  obs <- object$obs
  neq <- object$neq
  nvar <- object$nvar
  newnames <- colnames(newdata)
  is.rw <- !(winsize == 0)
  winsize <- abs (winsize)
  if(winsize > obs)
    winsize <- obs - 1
  totobs <- obs + n.ahead
  grid <- (1:totobs)/totobs
  object$z <- head(grid, obs)
  object$ez <- tail(grid, n.ahead)
  prediction <- matrix(NA, nrow = n.ahead, ncol = neq)
  colnames <- names(object$formula)
  X <- object$x
  Y <- object$y
  for (t in 1:n.ahead)
  {
    i <- ifelse (is.rw, obs - winsize, 1)
    object$y <- Y[i:(obs + t -1), ]
    object$z <- grid[i:(obs + t -1)]
    object$ez <- grid[obs + t]
    for (eq in 1:neq)
      object$x[[eq]] <- X[[eq]][i:(obs + t -1), ]
    theta <- try(update(object)$coefficients)
    if (inherits(theta, "try-error"))
    {
      object$bw <- bw (object)
      theta<- update(object)$coefficients
    }
    for (eq in 1:neq)
    {
      names <- colnames(object$x[[eq]])
      n.col <- nvar[eq]
      is.intercept <- "(Intercept)" %in% names
      n.col <- n.col - is.intercept
      index <- newnames%in%names
      if(sum(index) != n.col)
        stop("\nDataset 'newdata' does not contain all regressors in the model.\n")
      newx <- as.numeric(newdata[t, index])
      if(is.intercept)
        newx <- c(1L, newx) 
      columns <- (1 + ifelse(eq == 1, 0, sum(nvar[1:(eq-1)]))):sum(nvar[1:eq])
      prediction[t, eq] <- crossprod(theta[, columns], newx)
      X[[eq]] <- rbind(X[[eq]], newx)
    }
    Y <- rbind(Y, prediction[t, ])
  }
  colnames(prediction) <- names(object$formula)
  return(prediction)
}

#' @rdname forecast-tvReg
#' @method forecast tvplm
#' @examples
#' data(OECD)
#' tvpols <- tvPLM(lhe~lgdp+pop65+pop14+public, index = c("country", "year"), 
#' data = OECD, method = "pooling", bw =  8.9)
#' newdata <- OECD[c(7, 9), 4:7]
#' forecast(tvpols, newdata = newdata, n.ahead = 2)
#' @export
forecast.tvplm<-function (object, newdata, n.ahead = 1, winsize = 0, ...) 
{
  if(!inherits(object, c("tvplm")))
    stop("\nArgument 'object' should be entered and it should have class 'tvplm'.\n")
  if(!is.null(object$z))
    stop("\nYour model coefficients are functions of a random variable 'z', use function 
         'predict' with parameter 'newz'.\n")
  if(!inherits(newdata, c("pdata.frame", "matrix", "data.frame")))
    stop("\nArgument 'newdata' should be a 'pdata.frame', 
         'matrix' or 'data.frame'.\n")
  obs <- object$obs
  neq <- object$neq
  if(NROW(newdata) != n.ahead)
    stop("\nDimensions of 'newdata' are not compatible with 'n.ahead'.\n")
  if(NCOL(newdata) != object$nvar)
    stop("\nDimensions of 'newdata' are not compatible with number of variables in the model.\n")
  is.rw <- 1
  if(winsize != 0)
  {
    is.rw <- 2
    winsize <- abs (winsize)
    if(winsize > obs)
    {
      warning("\nArgument 'winsize' too big, set to abs(obs)-1: ", obs-1, "\n")
      winsize <- obs - 1
    }
  }
  else
    winsize <- obs
  method <- object$method
  prediction <- numeric(n.ahead)
  totobs <- obs + n.ahead
  grid <- (1:totobs)/totobs
  i <- obs - winsize + 1
  ind <- c(t(outer((0:(neq-1))*obs, i:obs, FUN = "+")))
  X <- object$x[ind,]
  Y <- object$y[ind]
  if(is.rw == 2)
    i <- obs - winsize 
  for (t in 1:n.ahead)
  {
    ind <- (i + (is.rw - 1) * t):(obs + t-1)
    object$z <- grid[ind]
    object$ez <- grid[obs + t]
    object$obs <- length(ind)
    object$x <- X
    object$y <- Y
    if(method != "within")
    {
      theta <- update(object)$coefficients
      if (inherits(theta, "try-error"))
      {
        object$bw <- bw(object)
        theta <- update(object)$coefficients
      }
    }
    else 
    {
      theta <- try(tvFE(object)$coefficients)
      if (inherits(theta, "try-error"))
      {
        object$bw <- bw(object)
        theta <- update(object)$coefficients
      }
    }
    regressor <- as.numeric(newdata[t,])
    prediction[t] <- theta%*%regressor
    X = Y <- NULL
    for (k in 1:neq)
    {
      ind <- (k-1)*object$obs + is.rw:object$obs
      X <- rbind(X, object$x[ind,], regressor)
      Y <- c(Y, object$y[ind], prediction[t])
    }
  }
  return(prediction)
}
  