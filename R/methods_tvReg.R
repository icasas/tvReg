
#' Coefficients and residuals of functions in tvReg
#' 
#' Return coefficients and residuals for objects with class attribute \code{tvlm}, 
#' \code{tvar}, \code{tvvar}, \code{tvirf}, \code{tvsure}.
#' 
#' @param object An object used to select a method.
#' @param ... Other parameters passed to specific methods.
#' @rdname methods.tvReg
#' @importFrom stats coef
#' @importFrom stats residuals
#' @export

coef.tvlm <- function(object, ...)
{
  return (object$tvcoef)
}
coef.tvvar <- function(object, ...)
{
  return (object$tvcoef)
}
coef.tvirf <- function(object, ...)
{
  return (object$tvcoef)
}
coef.tvsure <- function(object, ...)
{
  return (object$tvcoef)
}
residuals.tvlm <- function(object, ...)
{
  return (object$residuals)
}
residuals.tvvar <- function(object, ...)
{
  return (object$residuals)
}
residuals.tvirf <- function(object, ...)
{
  return (object$residuals)
}
residuals.tvsure <- function(object, ...)
{
  return (object$residuals)
}
