#' Forecast Methods for Objects in tvReg.
#'
#' Forecast methods for objects with class attribute \code{tvlm}, \code{tvar}, 
#' \code{tvvar}, \code{tvirf} and \code{tvsure}. If the 
#' smoothing variable (z) in the model is non-NULL and it is a random 
#' variable then use function \code{predict} with parameter \code{newz}.
#' @param object Object of class \code{tvlm}, \code{tvar}, \code{tvvar} or \code{tvsure}.
#' @param ... Other parameters passed to specific methods.
#' @return An object of class matrix or vector with the same dimensions than  the dependent 
#' variable of \code{object}.
#' @seealso \code{\link{predict}}.
#' @import methods
#' @export
forecast <- function(object, ...) UseMethod("forecast", object)

#' @rdname forecast
#' @method forecast tvlm
#' @param newx A vector, dataframe or matrix with new values of all variables in x. No need to 
#' input the intercept.
#' @param n.ahead A scalar with the forecast horizon, value 1 by default.
#' @param winsize A scalar. If 0 then an 'increase window' forecasting is performed.
#'  Otherwise a 'rolling window' forecasting is performed with window size given by 
#'  'winsize'.
#' @examples 
#' data("RV")
#' RV2 <- head(RV, 2001)
#' tvHAR <- tvLM (RV ~ RV_lag + RV_week + RV_month, data = RV2, bw = 20)
#' newx <- cbind(RV$RV_lag[2002:2004], RV$RV_week[2002:2004],
#'               RV$RV_month[2002:2004])
#' forecast(tvHAR, newx, n.ahead = 3)
#' 
#' @export
forecast.tvlm<-function (object, newx, n.ahead = 1, winsize = 0, ...) 
{
  if(!inherits(object, c("tvlm")))
    stop("\nParameter 'object' should be entered and it should have class 'tvlm' or 'tvar'.\n")
  if(!is.null(object$z))
    stop("\nYour model coefficients are functions of a random variable 'z', use function 
         'predict' with parameter 'newz'.\n")
  if(!inherits(newx, c("data.frame", "matrix", "numeric", "vector")))
    stop("\nParameter 'newx' should be entered and it should be a numeric vector if there is only
         one row or a 'matrix' or a 'data.frame' for more than one row.\n")
  if(n.ahead == 1) 
    newx <- matrix(newx, ncol = length(newx))
  if(NROW(newx) != n.ahead)
    stop("\nDimensions of 'newx' are not compatible with 'n.ahead'.\n")
  n.col <- NCOL(newx)
  is.intercept <- ("(Intercept)" %in% colnames(object$x))
  if(is.intercept & n.col ==(NCOL(object$x) - 1))
    newx <- cbind(rep(1, n.ahead), newx)
  obs <- object$obs
  is.rw <- !(winsize == 0)
  if(winsize > obs | winsize < 0)
    winsize <- obs - 1
  prediction <- numeric(n.ahead)
  tobs <- obs + n.ahead
  grid <- (1:tobs)/tobs
  X <- object$x
  Y <- object$y
  for (t in 1:n.ahead)
  {
    i <- ifelse (is.rw, obs - winsize, 1)
    object$x <- X[i:(obs + t -1),]
    object$y <- Y[i:(obs + t -1)]
    object$z <- grid[i:(obs + t -1)]
    object$ez <- grid[obs + t]
    beta <- try(tvOLS(object)$tvcoef)
    if (class (beta) == "try-error")
    {
      object$bw <- bw (object)
      beta<- tvOLS(object)$tvcoef
    }
    prediction[t] <- beta%*%newx[t,]
    X <- rbind(X, newx[t,])
    Y <- c (Y, prediction[t])
  }
  return(prediction)
}

#' @rdname forecast
#' @method forecast tvar
# @inheritParams forecast.tvlm
#' @param newz A vector with the new values of the smoothing variable.
#' @param newexogen A matrix or vector with the new values of the exogenous variables.
#' Only for predictions of *tvar* and *tvvar* objects.
#' @examples 
#' exogen = RV2[, c("RV_week", "RV_month")]
#' tvHAR2 <- tvAR(RV2$RV_lag, p = 1, exogen = exogen, bw = 20)
#' newexogen <- newx[, -1]
#' forecast(tvHAR2, n.ahead = 3, newexogen = newexogen)
#' 
#' @export
forecast.tvar <- function(object, n.ahead = 1, newz = NULL, newexogen = NULL, winsize = 0, ...) 
{
  if(!inherits(object, c("tvar")))
    stop("\nParameter 'object' should be entered and it should have class 'tvar'.\n")
  if(!inherits(newz, c( "numeric", "vector")) & !is.null(object$z))
    stop("\nArgument 'newz' should be entered and it should be a numeric vector.\n")
  if (!identical(NCOL(newexogen), NCOL(object$exogen))) 
    stop("\nWrong dimension in 'newexogen'.\n")
  is.exogen <- !is.null(newexogen)
  if(is.exogen)
  {
    if(NCOL(newexogen) == 1)
      newexogen <- matrix (newexogen, ncol = length(newexogen))
    if(NROW(newexogen) != n.ahead)
      stop("\nDimensions of 'newxexogen' and 'n.ahead' are not compatible.\n")
    newexogen <- as.matrix(newexogen)
  }
  is.intercept <- (object$type == "const")
  prediction <- numeric(n.ahead)
  obs <- object$obs
  is.rw <- !( winsize == 0)
  if(winsize > obs | winsize < 0)
    winsize <- obs - 1
  prediction <- numeric(n.ahead)
  tobs <- obs + n.ahead
  if(is.null(newz))
    grid <- (1:tobs)/tobs
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
    beta <- try(tvOLS(object)$tvcoef)
    if (class (beta) == "try-error")
    {
      object$bw <- bw (object)
      beta<- tvOLS(object)$tvcoef
    }
    newx <- tail(object$y, p)
    if(is.intercept) newx <- c(1L, newx)
    if(is.exogen) newx <- c(newx, newexogen[t, ])
    newx <- newx[mask]
    prediction[t] <- beta%*%newx
    X <- rbind(X, newx)
    Y <- c (Y, prediction[t])
  }
  return(prediction)
}


#' @rdname forecast
#' @method forecast tvvar
#' @inheritParams forecast.tvar
#' 
#' @examples 
#' data(usmacro, package = "bvarsv")
#' tvVAR <- tvVAR(usmacro, p = 6, type = "const", bw = c(1.8, 20, 20))
#' forecast(tvVAR, n.ahead = 10)
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
  if (!identical(NCOL(newexogen), NCOL(object$exogen))) 
    stop("\nWrong dimension in 'newexogen'.\n")
  is.exogen <- FALSE
  if(!is.null(newexogen))
  {
    is.exogen <- TRUE
    if(NCOL(newexogen) == 1)
      newexogen <- matrix (newexogen, ncol = length(newexogen))
    if(NROW(newexogen) != n.ahead)
      stop("\nDimensions of 'newxexogen' and 'n.ahead' are not compatible.\n")
    newexogen <- as.matrix(newexogen)
  }
  neq <- object$neq
  obs <- object$obs
  prediction <- matrix(NA, nrow = n.ahead, ncol = neq)
  colnames(prediction) <- colnames(object$y)
  is.rw <- !(winsize == 0)
  if(winsize > obs | winsize < 0)
    winsize <- obs - 1
  tobs <- obs + n.ahead
  if(is.null(newz))
    grid <- (1:tobs)/tobs
  else
    grid <- c(object$z, newz)
  X <- object$x
  Y <- object$y
  p <- object$p
  nlags <- p*neq
  for (t in 1:n.ahead)
  {
    i <- ifelse (is.rw, obs - winsize, 1)
    object$x <- X[i:(obs + t -1), ]
    object$y <- Y[i:(obs + t -1), ]
    object$z <- grid[i:(obs + t -1)]
    object$ez <- grid[obs + t]
    beta <- try(tvOLS(object)$tvcoef)
    if (class (beta) == "try-error")
    {
      object$bw <- bw (object)
      beta<- tvOLS(object)$tvcoef
    }
    rhs <- c(tail(object$y, 1))
    rhs <- c(rhs, tail(object$x, 1)[1:(nlags - neq)])
    if(is.intercept) 
      rhs <- c(rhs, "(Intercept)" = 1L)
    if(is.exogen)
      rhs <- c(rhs, newexogen[t, ])
    for (i in 1:neq)
      prediction[t, i] <- beta[[i]]%*%rhs
    X <- rbind(X, rhs)
    Y <- rbind(Y, prediction[t, ])
  }
  return(prediction)
}


#' @rdname forecast
#' @method forecast tvsure
# @inheritParams forecast.tvlm
#' @param newdata A matrix or data.frame with the values of the regressors to use
#' for forecasting.
#' @examples 
#' data("Kmenta", package = "systemfit")
#' eqDemand <- consump ~ price + income
#' eqSupply <- consump ~ price + farmPrice 
#' system <- list(demand = eqDemand, supply = eqSupply)
#' tvOLS.fit <- tvSURE(system, data = Kmenta, est = "ll", bw = c(1.5, 1.5))
#' newdata <- data.frame(consump = c(95, 100, 102), price = c(90, 100, 103), 
#' farmPrice = c(70, 95, 103), income = c(82, 94, 115))
#' forecast(tvOLS.fit, newdata = newdata, n.ahead = 3)
#' 
#' @export
forecast.tvsure<-function (object, newdata, n.ahead = 1, winsize = 0, ...) 
{
  if(!inherits(object, c("tvsure")))
    stop("\nParameter 'object' should be entered and it should have class 'tvsure'.\n")
  if(!is.null(object$z))
    stop("\nYour model coefficients are functions of a random variable 'z', use function 
         'predict' with parameter 'newz'.\n")
  if(!inherits(newdata, c("matrix", "data.frame")))
    stop("\nParameter 'newdata' should be entered and it should be a 'matrix' or a 'data.frame'.\n")
  newdata <- as.matrix(newdata)
  if(NROW(newdata) != n.ahead)
    stop("\nDimensions of 'newdata' are not compatible with 'n.ahead'.\n")
  obs <- object$obs
  neq <- object$neq
  nvar <- object$nvar
  newnames <- colnames(newdata)
  is.rw <- !(winsize == 0)
  if(winsize > obs | winsize < 0)
    winsize <- obs - 1
  tobs <- obs + n.ahead
  grid <- (1:tobs)/tobs
  X <- object$x
  Y <- object$y
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
    beta <- try(update(object)$tvcoef)
    if (class (beta) == "try-error")
    {
      object$bw <- bw (object)
      beta<- update(object)$tvcoef
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
      prediction[t, eq] <- beta[, columns]%*%newx
      X[[eq]] <- rbind(X[[eq]], newx)
    }
    Y <- rbind(Y, prediction[t, ])
  }
  colnames(prediction) <- names(object$formula)
  return(prediction)
}
  