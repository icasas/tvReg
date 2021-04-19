
#' Predict Methods for Objects in tvReg.
#'
#' Predict methods for objects with class attribute \code{tvlm}, \code{tvar}, 
#' \code{tvvar}, \code{tvirf}, \code{tvsure} and \code{tvplm}. This function needs new values of 
#' variables y (response), x (regressors), exogen (exogenous variables, when used),
#' and  z (smoothing variable). 
#' @param object An object used to select a method.
#' @param newdata A data.frame, matrix or vector with new values of all independent variables. No need to 
#' input the intercept. Note that this does not refer to the variables in 'exogen' which might be part of the 'tvar' 
#' and 'tvvar' objects. Those must be included in 'newexogen'.
#' @param newdata A vector or matrix with new values of the lags included in the model.
#' @param newz A vector with new values of the smoothing variable.
#' @param newexogen A matrix or vector with the new value of the exogenous variables.
#' Only for predictions of 'tvar' and 'tvvar' objects.
#' @param ... Other arguments passed to specific methods.
#' @return An object of class matrix or vector with the prediction.
#' @rdname predict-tvReg
#' @method predict tvlm
#' @seealso \code{\link{forecast}}.
#' @examples
#' ## Example of TVLM prediction with coefficients as 
#' ## functions of the realized quarticity
#' 
#' data("RV")
#' RV2 <- head(RV, 2001)
#' z <- RV2$RQ_lag_sqrt
#' TVHARQ <- tvLM (RV ~ RV_lag + RV_week + RV_month, 
#'                  z = z, data = RV2, bw = 0.0062)
#' newdata <- cbind(RV$RV_lag[2002:2004], RV$RV_week[2002:2004],
#'               RV$RV_month[2002:2004])
#' newz <- RV$RQ_lag_sqrt[2002:2004]
#' predict(TVHARQ, newdata, newz)
#' 
#' @export
predict.tvlm<-function (object, newdata, newz, ...) 
{
  if(!inherits(object, c("tvlm")))
    stop("\nArgument 'object' should be entered and it should have class 'tvlm'.\n")
  if(is.null(newdata))
    return(stats::fitted(object, ...))
  if(is.null(object$z))
    stop("\nYour model coefficients are functions of time, use function 
         'forecast' with argument 'n.ahead' as horizon.\n")
  if(!inherits(newdata, c("data.frame", "matrix", "numeric", "vector")))
    stop("\nArgument 'newdata' should a numeric vector if there is only
         one row or a 'matrix' or a 'data.frame' for more than one row.\n")
  if(!inherits(newz, c( "numeric", "vector")))
    stop("\nArgument 'newz' should be entered and it should be a numeric vector.\n")
  n.ahead <- length(newz)
  if(n.ahead == 1) 
    newdata <- matrix(newdata, ncol = length(newdata))
  n.col <- NCOL(newdata)
  if(!identical(NROW(newdata), n.ahead))
    stop("\nDimensions of 'newdata' and 'newz' are not compatible\n")
  obs <- NROW(object$x) 
  is.intercept <- ("(Intercept)" %in% colnames(object$x))
  if(is.intercept & n.col == (NCOL(object$x) - 1))
    newdata <- cbind(rep(1, n.ahead), newdata)
  prediction <- numeric(n.ahead)
  object$ez <- newz
  theta <- tvOLS(object)$coefficients
  for (t in 1:n.ahead)
    prediction[t] <- theta[t, ]%*%newdata[t, ]
  return(prediction)
}

#' @rdname predict-tvReg
#' @method predict tvar
#' @examples 
#' ## Example of TVAR prediction with coefficients as 
#' ## functions of the realized quarticity
#' 
#' exogen = RV2[, c("RV_week", "RV_month")]
#' TVHARQ2 <- tvAR (RV2$RV, p = 1, exogen = exogen,  
#'                       z = RV2[, "RQ_lag_sqrt"], bw = 0.0062)
#' newylag <- RV$RV[2002:2004]
#' newz <- RV$RQ_lag_sqrt[2002:2004]
#' newexogen <- RV[2002:2004, c("RV_week", "RV_month")]
#' predict(TVHARQ2, newylag,  newz, newexogen = newexogen)
#' @export
predict.tvar<-function (object, newdata, newz, newexogen = NULL, ...) 
{
  if(!inherits(object, c("tvar")))
    stop("\nArgument 'object' should be entered and it should have class 'tvar'.\n")
  if(is.null(newdata))
    return(stats::fitted(object, ...))
  if(is.null(object$z))
    stop("\nYour model coefficients are functions of time, use function 
         'forecast' with argument 'n.ahead' as horizon.\n")
  if(!inherits(newdata, c("vector", "numeric", "data.frame", "matrix")))
    stop("\nArgument 'newdata' should be a numeric vector or matrix.\n")
  if(!inherits(newz, c( "numeric", "vector")))
    stop("\nArgument 'newz' should be entered and it should be a numeric vector.\n")
  if (!identical(NCOL(newexogen), NCOL(object$exogen))) 
    stop("\nWrong dimension in 'newexogen'.\n")
  if(!identical(NROW(newdata), NROW(newz)))
    stop("\nDifferent number or row in 'newdata' and 'newz'.\n")
  if(!identical(NROW(newdata), NROW(newexogen)) & !is.null(newexogen))
    stop("\nDifferent number or row in 'newdata' and 'newexogen'.\n")
  n.ahead <- length(newz)
  newx <- newdata
  p <- object$p
  if(p == 1)
    newx <- matrix(newx, ncol = 1)
  mask <- object$mask
  if(!is.null(newexogen))
  {
    if(NCOL(newexogen) == 1)
      newexogen <- matrix (newexogen, ncol = length(newexogen))
    if(NROW(newexogen) != n.ahead)
      stop("\nDimensions of 'newxexogen' and 'n.ahead' are not compatible.\n")
    newexogen <- as.matrix(newexogen)
    newx <- cbind(newx, newexogen)
  }
  n.col <- NCOL(newx)
  if(!identical(NROW(newx), n.ahead))
    stop("\nDimensions of 'newdata' and 'newz' are not compatible\n")
  obs <- NROW(object$x) 
  is.intercept <- ("(Intercept)" %in% colnames(object$x))
  if(is.intercept & n.col == (NCOL(object$x) - 1))
    newx <- cbind(rep(1, n.ahead), newx)
  prediction <- numeric(n.ahead)
  object$ez <- newz
  theta <- tvOLS(object)$coefficients
  for (t in 1:n.ahead)
    prediction[t] <- theta[t, ]%*%newx[t, mask]
  return(prediction)
}
#' @rdname predict-tvReg
#' @method predict tvvar
#' @examples 
#' ## Example of TVVAR prediction with coefficients as 
#' ## functions of a random ARMA (2,2) process
#' 
#' data(usmacro, package = "bvarsv")
#' smoothing <- arima.sim(n = NROW(usmacro) + 3, 
#' list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)), 
#' sd = sqrt(0.1796))
#' smoothing <- as.numeric(smoothing)
#' TVVAR.z <- tvVAR(usmacro, p = 6, type = "const", 
#'                z = smoothing[1:NROW(usmacro)], bw = c(16.3, 16.3, 16.3))
#' newdata <- data.frame(inf = c(2, 1, 6), une = c(5, 4, 9), tbi = c(1, 2.5, 3))
#' newz <- c(0, 1.2, -0.2)
#' predict(TVVAR.z, newdata = newdata, newz = newz)
#' 
#' @export
predict.tvvar<-function (object, newdata, newz, newexogen = NULL, ...) 
{
  if(!inherits(object, c("tvvar")))
    stop("\nArgument 'object' should be entered and it should have class 'tvvar'.\n")
  if(is.null(newdata))
    return(stats::fitted(object, ...))
  if(is.null(object$z))
    stop("\nYour model coefficients are functions of time, use function 
         'forecast' with argument 'n.ahead' as horizon.\n")
  if(!inherits(newdata, c("vector", "numeric", "data.frame", "matrix")))
    stop("\nArgument 'newdata' should be a numeric vector if there is only
         one row or a 'matrix' or a 'data.frame' for more than one row.\n")
  if(!inherits(newz, c( "numeric", "vector")))
    stop("\nArgument 'newz' should be entered and it should be a numeric vector.\n")
  if(!identical(NROW(newdata), NROW(newz)))
    stop("\nDifferent number or row in 'newdata' and 'newz'.\n")
  if(!identical(NROW(newdata), NROW(newexogen)) & !is.null(newexogen))
    stop("\nDifferent number or row in 'newdata' and 'newexogen'.\n")
  n.ahead <- length(newz)
  neq <- object$neq
  p <- object$p
  newx <- newdata
  nlags <- neq*p
  if(neq == 1)
    newx <- matrix(newdata, ncol = 1)
  newx <- as.matrix(newx)
  if(NROW(newx) != n.ahead)
    stop("\nDimensions of 'newx' and 'newz' are not compatible\n")
  if (!identical(NCOL(newexogen), NCOL(object$exogen))) 
    stop("\nWrong dimension in 'newexogen'.\n")
  is.exogen <- FALSE
  is.intercept <- ifelse(object$type == "const", TRUE, FALSE)
  if(!is.null(newexogen))
  {
    is.exogen <- TRUE
    if(NCOL(newexogen) == 1)
      newexogen <- matrix (newexogen, ncol = length(newexogen))
    if(NROW(newexogen) != n.ahead)
      stop("\nDimensions of 'newxexogen' and 'n.ahead' are not compatible.\n")
    n.exogen <- NCOL(newexogen)
    newexogen <- as.matrix(newexogen)
  }
  prediction <- matrix(NA, nrow = n.ahead, ncol = neq)
  colnames(prediction) <- colnames(object$y)
  object$ez <- newz
  theta <- tvOLS(object)$coefficients
  rhs <- tail(object$x, 1)
  rhs <- as.numeric(rhs)
  for (t in 1:n.ahead)
  {
    rhs <- c(newx[t, ], rhs[1:(nlags - neq)])
    if(is.intercept) 
      rhs <- c(rhs, "(Intercept)" = 1L)
    if(is.exogen)
      rhs <- c(rhs, newexogen[t, ])
    for (i in 1:neq)
      prediction[t, ] <- theta[[i]][t,]%*%rhs
  }
  return(prediction)
}

#' @rdname predict-tvReg
#' @method predict tvsure
#' @param newdata A dataframe with new values of all regressors, with the
#' same name and order as they appear in argument 'data' from the 'tvsure'
#' object
#' @examples 
#' ## Example of TVSURE prediction with coefficients as 
#' ## functions of an ARMA(2,2) process
#' data("Kmenta", package = "systemfit")
#' nobs <- NROW (Kmenta)
#' eqDemand <- consump ~ price + income
#' eqSupply <- consump ~ price + farmPrice 
#' system <- list(demand = eqDemand, supply = eqSupply)
#' smoothing <- arima.sim(n = nobs + 3, 
#'                        list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)), 
#'                        sd = sqrt(0.1796))
#' smoothing <- as.numeric(smoothing)
#' TVOLS.z <- tvSURE(system, data = Kmenta,  
#'                       z = smoothing[1:nobs],  bw = c(7, 1.8),
#'                       est = "ll")
#' newdata <- data.frame(consump = c(95, 100, 102), price = c(90, 100, 103), 
#'                       farmPrice = c(70, 95, 103), income = c(82, 94, 115))
#' newz <- tail(smoothing, 3)
#' predict(TVOLS.z, newdata = newdata, newz = newz)
#' 
#' @export
predict.tvsure<-function (object, newdata, newz, ...) 
{
  if(!inherits(object, c("tvsure")))
    stop("\nArgument 'object' should be entered and it should have class 'tvsure'.\n")
  if(is.null(newdata))
    return(stats::fitted(object, ...))
  if(is.null(object$z))
    stop("\nYour model coefficients are functions of time, use function 
         'forecast' with argument 'n.ahead' as horizon.\n")
  if(!inherits(newdata, c("matrix", "data.frame")))
    stop("\nArgument 'newdata' should be a 'matrix' or a 'data.frame'.\n")
  if(!inherits(newz, c( "numeric", "vector")))
    stop("\nArgument 'newz' should be entered and it should be a numeric vector.\n")
  newdata <- as.matrix(newdata)
  neq <- object$neq
  nvar <- object$nvar
  n.ahead <- length(newz)
  if(NROW(newdata) != n.ahead)
    stop("\nDimensions of 'newdata' and 'newz' are incompatible.\n")
  newx <- list()
  newnames <- colnames(newdata)
  object$ez <- newz
  theta <- update(object)$coefficients
  prediction <- matrix(NA, nrow = n.ahead, ncol = neq)
  for (t in 1:n.ahead)
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
      prediction[t, eq] <- theta[t, columns]%*%newx
    }
  colnames(prediction) <- names(object$formula)
  return(prediction)
}

#' @rdname predict-tvReg
#' @method predict tvplm
#' @param newdata A pdata.frame with new values of all regressors, with the
#' same name and order as they appear in argument 'data' from the 'tvplm'
#' object
#' @examples
#' data(OECD)
#' z <- runif(length(levels(OECD$year)), 10, 15)
#' TVPOLS <- tvPLM(lhe~lgdp+pop65+pop14+public, z = z,
#' index = c("country", "year"), data = OECD,  method ="pooling", bw =  2)
#' newdata <- cbind(lgdp = c(10, 13), pop65 = c(9, 12), 
#' pop14 = c(17, 30), public = c(13, 20))  
#' newz <- runif(2, 10, 15)
#' predict(TVPOLS, newdata = newdata, newz = newz)
#' @export
predict.tvplm<-function (object, newdata, newz, ...) 
{
  if(!inherits(object, c("tvplm")))
    stop("\nArgument 'object' should be entered and it should have class 'tvplm'.\n")
  if(is.null(newdata))
    return(stats::fitted(object, ...))
  if(is.null(object$z))
    stop("\nYour model coefficients are functions of time, use function 
         'forecast' with argument 'n.ahead' as horizon.\n")
  if(!inherits(newdata, c("pdata.frame", "matrix", "data.frame")))
    stop("\nArgument 'newdata' should be a 'pdata.frame', 
         'matrix' or 'data.frame'.\n")
  if(!inherits(newz, c( "numeric", "vector")))
    stop("\nArgument 'newz' should be a numeric vector.\n")
  n.ahead <- length(newz)
  if(NROW(newdata) != n.ahead)
    stop("\nDimensions of 'newdata' and 'newz' are incompatible.\n")
  obs <- object$obs
  neq <- object$neq
  object$ez <- newz
  object.up <- update(object)
  theta <- object.up$coefficients
  method <- object$method
  if(method != "within")
  {
    prediction <- numeric(n.ahead)
    for (t in 1:n.ahead)
      prediction[t] <- crossprod(theta[t,], newdata[t,])
    return(prediction)
  }
  prediction <- matrix(NA, nrow = n.ahead, ncol = neq)
  for (t in 1:n.ahead)
    prediction[t, ] <- crossprod(theta[t,], newdata[t,])
  prediction <- sweep(prediction, 2, object.up$alpha, "+")
  colnames(prediction) <- levels (object$index[, 1])
  #Fitted and residuals must be calculated here for bootstrap
  return(prediction)
}