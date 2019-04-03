
#' Predict Methods for Objects in tvReg.
#'
#' Predict methods for objects with class attribute \code{tvlm}, \code{tvar}, 
#' \code{tvvar}, \code{tvirf} and \code{tvsure}. This function needs new values of 
#' variables y (response), x (regressors), exogen (exogenous variables, when used),
#' and  z (smoothing variable). 
#' @param object Object of class \code{tvlm}, \code{tvar}, \code{tvvar} or \code{tvsure}.
#' @param newy A vector with new values of the response variable 
#' @param newx A dataframe with new values of all variables in x. No need to 
#' input the intercept.
#' @param newz A vector with new values of the smoothing variable.
#' @param newexogen A matrix or vector with the new value of the exogenous variables.
#' Only for predictions of 'tvar' and 'tvvar' objects.
#' @param ... Other arguments passed to specific methods.
#' @return An object of class matrix or vector with the prediction.
#' @rdname predict-tvReg
#' @method predict tvlm
#' @seealso \code{\link{forecast}}.
#' @examples
#' ## Example of TV-LM prediction with coefficients as 
#' ## functions of the realized quarticity
#' 
#' data("RV")
#' RV2 <- head(RV, 2001)
#' z <- RV2$RQ_lag_sqrt
#' tvHARQ <- tvLM (RV ~ RV_lag + RV_week + RV_month, 
#'                  z = z, data = RV2, bw = 0.0062)
#' newx <- cbind(RV$RV_lag[2002:2004], RV$RV_week[2002:2004],
#'               RV$RV_month[2002:2004])
#' newz <- RV$RQ_lag_sqrt[2002:2004]
#' predict(tvHARQ, newx, newz)
#' 
#' @export
predict.tvlm<-function (object, newx, newz, ...) 
{
  if(!inherits(object, c("tvlm")))
    stop("\nArgument 'object' should be entered and it should have class 'tvlm'.\n")
  if(is.null(object$z))
    stop("\nYour model coefficients are functions of time, use function 
         'forecast' with argument 'n.ahead' as horizon.\n")
  if(!inherits(newx, c("data.frame", "matrix", "numeric", "vector")))
    stop("\nArgument 'newx' should be entered and it should be a numeric vector if there is only
         one row or a 'matrix' or a 'data.frame' for more than one row.\n")
  if(!inherits(newz, c( "numeric", "vector")))
    stop("\nArgument 'newz' should be entered and it should be a numeric vector.\n")
  n.ahead <- length(newz)
  if(n.ahead == 1) 
    newx <- matrix(newx, ncol = length(newx))
  n.col <- NCOL(newx)
  if(!identical(NROW(newx), n.ahead))
    stop("\nDimensions of 'newx' and 'newz' are not compatible\n")
  obs <- NROW(object$x) 
  is.intercept <- ("(Intercept)" %in% colnames(object$x))
  if(is.intercept & n.col == (NCOL(object$x) - 1))
    newx <- cbind(rep(1, n.ahead), newx)
  prediction <- numeric(n.ahead)
  object$ez <- newz
  beta <- tvOLS(object)$tvcoef
  for (t in 1:n.ahead)
    prediction[t] <- beta[t, ]%*%newx[t, ]
  return(prediction)
}

#' @rdname predict-tvReg
#' @method predict tvar
#' @examples 
#' ## Example of TV-AR prediction with coefficients as 
#' ## functions of the realized quarticity
#' 
#' exogen = RV2[, c("RV_week", "RV_month")]
#' tvHARQ2 <- tvAR (RV2$RV, p = 1, exogen = exogen,  
#'                       z = RV2[, "RQ_lag_sqrt"], bw = 0.0062)
#' newylag <- RV$RV[2002:2004]
#' newz <- RV$RQ_lag_sqrt[2002:2004]
#' newexogen <- RV[2002:2004, c("RV_week", "RV_month")]
#' predict(tvHARQ2, newylag,  newz, newexogen = newexogen)
#' @export
predict.tvar<-function (object, newy, newz, newexogen = NULL, ...) 
{
  if(!inherits(object, c("tvar")))
    stop("\nArgument 'object' should be entered and it should have class 'tvar'.\n")
  if(is.null(object$z))
    stop("\nYour model coefficients are functions of time, use function 
         'forecast' with argument 'n.ahead' as horizon.\n")
  if(!inherits(newy, c("numeric", "vector")))
    stop("\nArgument 'newy' should be entered and it should be a numeric vector.\n")
  if(!inherits(newz, c( "numeric", "vector")))
    stop("\nArgument 'newz' should be entered and it should be a numeric vector.\n")
  if (!identical(NCOL(newexogen), NCOL(object$exogen))) 
    stop("\nWrong dimension in 'newexogen'.\n")
  n.ahead <- length(newz)
  newx <- newy
  p <- object$p
  if(p == 1)
    newx <- matrix(newx, ncol = 1)
  newx <- as.matrix(newx)
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
    stop("\nDimensions of 'newx' and 'newz' are not compatible\n")
  obs <- NROW(object$x) 
  is.intercept <- ("(Intercept)" %in% colnames(object$x))
  if(is.intercept & n.col == (NCOL(object$x) - 1))
    newx <- cbind(rep(1, n.ahead), newx)
  prediction <- numeric(n.ahead)
  object$ez <- newz
  beta <- tvOLS(object)$tvcoef
  for (t in 1:n.ahead)
    prediction[t] <- beta[t, ]%*%newx[t, ]
  return(prediction)
}
#' @rdname predict-tvReg
#' @method predict tvvar
#' @examples 
#' ## Example of TV-VAR prediction with coefficients as 
#' ## functions of a random ARMA (2,2) process
#' 
#' data(usmacro, package = "bvarsv")
#' smoothing <- arima.sim(n = nrow(usmacro) + 3, 
#' list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)), 
#' sd = sqrt(0.1796))
#' smoothing <- as.numeric(smoothing)
#' tvVAR.z <- tvVAR(usmacro, p = 6, type = "const", 
#'                z = smoothing[1:nrow(usmacro)], bw = c(16.3, 16.3, 16.3))
#' newy <- data.frame(inf = c(2, 1, 6), une = c(5, 4, 9), tbi = c(1, 2.5, 3))
#' newz <- c(0, 1.2, -0.2)
#' predict(tvVAR.z, newy = newy, newz = newz)
#' 
#' @export
predict.tvvar<-function (object, newy, newz, newexogen = NULL, ...) 
{
  if(!inherits(object, c("tvvar")))
    stop("\nArgument 'object' should be entered and it should have class 'tvvar'.\n")
  if(is.null(object$z))
    stop("\nYour model coefficients are functions of time, use function 
         'forecast' with argument 'n.ahead' as horizon.\n")
  if(!inherits(newy, c("vector", "numeric", "data.frame", "matrix")))
    stop("\nArgument 'newy' should be entered and it should be a numeric vector if there is only
         one row or a 'matrix' or a 'data.frame' for more than one row.\n")
  if(!inherits(newz, c( "numeric", "vector")))
    stop("\nArgument 'newz' should be entered and it should be a numeric vector.\n")
  n.ahead <- length(newz)
  neq <- object$neq
  p <- object$p
  newx <- newy
  nlags <- neq*p
  is.exogen = FALSE
  is.intercept = ifelse(object$type == "const", TRUE, FALSE)
  if(neq == 1)
    newx <- matrix(newy, ncol = 1)
  newx <- as.matrix(newx)
  if(NROW(newx) != n.ahead)
    stop("\nDimensions of 'newx' and 'newz' are not compatible\n")
  if (!identical(NCOL(newexogen), NCOL(object$exogen))) 
    stop("\nWrong dimension in 'newexogen'.\n")
  if(!is.null(newexogen))
  {
    is.exogen = TRUE
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
  beta <- tvOLS(object)$tvcoef
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
      prediction[t, ] <- beta[[i]][t,]%*%rhs
  }
  return(prediction)
}

#' @rdname predict-tvReg
#' @method predict tvsure
#' @param newdata A dataframe with new values of all regressors, with the
#' same name and order as they appear in argument 'data' from the 'tvsure'
#' object
#' @examples 
#' ## Example of TV-SURE prediction with coefficients as 
#' ## functions of an ARMA(2,2) process
#' data("Kmenta", package = "systemfit")
#' nobs <- nrow (Kmenta)
#' eqDemand <- consump ~ price + income
#' eqSupply <- consump ~ price + farmPrice 
#' system <- list(demand = eqDemand, supply = eqSupply)
#' smoothing <- arima.sim(n = nobs + 3, 
#'                        list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)), 
#'                        sd = sqrt(0.1796))
#' smoothing <- as.numeric(smoothing)
#' tvOLS.z.fit <- tvSURE(system, data = Kmenta,  
#'                       z = smoothing[1:nobs],  bw = c(7, 1.8),
#'                       est = "ll")
#' newdata <- data.frame(consump = c(95, 100, 102), price = c(90, 100, 103), 
#'                       farmPrice = c(70, 95, 103), income = c(82, 94, 115))
#' newz <- tail(smoothing, 3)
#' predict(tvOLS.z.fit, newdata = newdata, newz = newz)
#' 
#' @export
predict.tvsure<-function (object, newdata, newz, ...) 
{
  if(!inherits(object, c("tvsure")))
    stop("\nArgument 'object' should be entered and it should have class 'tvsure'.\n")
  if(is.null(object$z))
    stop("\nYour model coefficients are functions of time, use function 
         'forecast' with argument 'n.ahead' as horizon.\n")
  if(!inherits(newdata, c("matrix", "data.frame")))
    stop("\nArgument 'newdata' should be entered and it should be a 'matrix' or a 'data.frame'.\n")
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
  beta <- update(object)$tvcoef
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
      prediction[t, eq] <- beta[t, columns]%*%newx
    }
  colnames(prediction) <- names(object$formula)
  return(prediction)
}