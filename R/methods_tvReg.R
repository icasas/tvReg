
#' @name tvReg-methods
#' @aliases coefficients.tvlm, coefficients.tvar,  coefficients.tvvar,  
#' coefficients.tvirf,  coefficients.tvsure
#' @title Coefficients and residuals of functions in tvReg
#' @description Return coefficients and residuals for objects with class attribute \code{tvlm}, 
#' \code{tvar}, \code{tvvar}, \code{tvirf}, \code{tvsure}.
#' @param object An object used to select a method.
#' @param ... Other parameters passed to specific methods.
#' @importFrom stats coef
#' 
#' @method coef tvlm
#' @export
coef.tvlm <- function(object, ...) return(object$tvcoef)

#' @rdname tvReg-methods
#' @method coef tvar
#' @export
coef.tvar <- function(object, ...) return(object$tvcoef)

#' @rdname tvReg-methods
#' @method coef tvvar
#' @export
coef.tvvar <- function(object, ...) return(object$tvcoef)

#' @rdname tvReg-methods
#' @method coef tvirf
#' @export
coef.tvirf <- function(object, ...) return(object$tvcoef)

#' @rdname tvReg-methods
#' @method coef tvsure
#' @export
coef.tvsure <- function(object, ...) return(object$tvcoef)

