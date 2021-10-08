#' Update and Re-fit the Models of package tvReg
#'
#' Update and Re-fit the Models of package tvReg
#' @rdname update.tvReg
#' @param object An object of any class in package tvReg.
#' @param ... Other parameters passed to specific methods.
#' @return An object of the same class than the argument *object*.
#' @method update tvlm
#' @export
update.tvlm <- function(object, ...)
{
  result <- tvOLS(object)
  object$coefficients <- result$coefficients
  object$fitted <- result$fitted
  object$residuals <- result$residuals
  return(object)
}

#' @rdname update.tvReg
#' @method update tvar 
#' @export 
update.tvar <- update.tvlm


#' @rdname update.tvReg
#' @method update tvvar
#' @export
update.tvvar <- function(object, ...)
{
  results <- tvOLS(object)
  object$coefficients <- results$coefficients
  object$fitted <- results$fitted
  object$residuals <- results$residuals
  class(object) <- "tvvar"
  return(object)
}
#' @rdname update.tvReg
#' @method update tvsure
#' @export
update.tvsure <- function(object, ...)
{
  method <- object$method
  tkernel <- object$tkernel
  obs <- object$obs
  neq <- object$neq
  if (method %in% c("identity", "tvOLS"))
    object$Sigma <- NULL
  result <- tvGLS(object)
  object$coefficients <- result$coefficients
  class(object) <- "tvsure"
  return(object)
}

#' @rdname update.tvReg
#' @method update tvplm
#' @export
update.tvplm <- function(object, ...)
{
  method <- object$method
  if (method != "within")
  {
    object$Sigma <- NULL
    result <- tvRE(object)
    if(method == "random")
    {
      result <- tvRE(object)
      itertemp <- 1
      tol <- object$control$tol
      maxiter <- object$control$maxiter
      tolold <- sum(result$coefficients^2)
      tolnew <- 0
      while((abs(tolold-tolnew)>tol) && (itertemp < maxiter))
      {
        tolold <- tolnew
        object$Sigma <- result$Sigma
        temp <- tvRE(object)
        tolnew <- sqrt(sum((result$coefficients - temp$coefficients)^2)/sum(result$coefficients^2))
        result <- temp
        itertemp <- itertemp + 1
      }
    }
  }
  else 
    result <- tvFE(object)
  object$coefficients <- result$coefficients
  object$alpha <- result$alpha
  object$fitted <- result$fitted
  object$residuals <- result$residuals
  object$alpha <- result$alpha
  class(object) <- "tvplm"
  return(object)
}