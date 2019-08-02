#' Update and Re-fit the Models of package tvReg
#'
#' Update and Re-fit the Models of package tvReg
#' @rdname update.tvReg
#' @param object An object of any class in package tvReg
#' @param ... Other parameters passed to specific methods.
#' @return An object of the same class than the argument *object*.
#' @method update tvlm
#' @export
update.tvlm <- function(object, ...)
{
  result <- tvOLS(object)
  object$tvcoef <- result$tvcoef
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
  object$tvcoef <- results$tvcoef
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
  if (method %in% c("identity", "tvOLS", "tvFGLS"))
  {
    object$Sigma <- NULL
    result <- tvGLS(object)
    Sigma <- array(rep(crossprod(result$residuals)/(obs - neq), obs), dim = c(neq, neq, obs))
    object$Sigma <- Sigma
  }
  if(method == "tvFGLS")
  {
    bw.cov <- object$bw.cov
    Sigma <- tvCov(x = result$residuals, bw = bw.cov, tkernel = tkernel)
    object$Sigma <- Sigma
    result <- tvGLS(object)
    itertemp <- 0
    tol <- object$control$tol
    maxiter <- object$control$maxiter
    tolold <- sum(result$tvcoef^2)
    tolnew <- 0
    while((abs(tolold-tolnew)>tol) && (itertemp <= maxiter))
    {
      tolold <- tolnew
      Sigma <- tvCov(x = result$residuals, bw = bw.cov, tkernel = tkernel)
      object$Sigma <- Sigma
      temp <- tvGLS(object)
      tolnew <- sqrt(sum((result$tvcoef - temp$tvcoef)^2)/sum(result$tvcoef^2))
      result <- temp
      itertemp <- itertemp + 1
    }
  }
  else if(method == "tvGLS")
  {
    result <- tvGLS(object)
  }
  object$tvcoef <- result$tvcoef
  object$fitted <- result$fitted
  object$residuals <- result$residuals
  class(object) <- "tvsure"
  return(object)
}