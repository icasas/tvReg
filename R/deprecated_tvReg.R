#' Deprecated function(s) in the tvReg package
#' 
#' These functions are provided for compatibility with older version of
#' the tvReg package.  They may eventually be completely removed.
#' @rdname tvReg-deprecated
#' @param object Object of class \code{tvsure}, class \code{tvvar} or class \code{tvirf}.
#' @param level Numeric, the confidence level required (between 0 and 1). 
#' @param runs (optional) Number of bootstrap replications.
#' @param tboot Type of wild bootstrap, choices 'wild'(default), 'wild2'. Option 'wild' uses the
#' distribution suggested by Mammen (1993) in the wild resampling, while 'wild2' uses the standard
#' normal.
#' @param ... Other parameters passed to specific methods.
#
#' @aliases CI confint
#' @section Details:
#' \tabular{rl}{
#'   \code{CI} \tab now a synonym for \code{\link{confint}}
#' }
#' 
#' @export
CI<- function(object, ...) 
{
  .Deprecated("confint", package = "tvReg")
  UseMethod("confint", object)
}

#' @rdname tvReg-deprecated
#' @method CI default
#' @export 
CI.default <- function(object, level = 0, runs = 0, tboot = NULL, ...)
{
  .Deprecated("confint", package = "tvReg")
  confint (object, level = level, runs = runs, tboot = tboot)
}

#' @rdname tvReg-deprecated
#' @method CI tvirf
#' @export
CI.tvirf <- function(object, level = 0, runs = 0, tboot = NULL , ...)
{
  .Deprecated("confint.tvirf", package = "tvReg")
  confint (object, level = level, runs = runs, tboot = tboot)
}
