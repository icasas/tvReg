##' Update and Re-fit the Models of package tvReg
##'
##' Update and Re-fit the Models of package tvReg
##' @rdname update.tvReg
##' @param object An object used to select a method.
##' @param ... Other parameters passed to specific methods.
##' @param y The dependent variable to update the model.
##' @return An object of class \code{tvsure}.
##' @method update tvsure
##' @export
update.tvsure <- function(object, y = NULL, ...)
{
  neq <- object$neq
  obs <- object$obs
  method <- object$method
  if(!is.null(y))
    object$y <- y
  if (method == "tvOLS" | method == "identity")
  {
    Sigma.ols <- array(rep(diag(1, neq), obs), dim=c(neq,neq, obs))
    object$Sigma <- Sigma.ols
  }
  result <- tvGLS(object)
  object$tvcoef <- result$tvcoef
  object$fitted <- result$fitted
  object$residuals <- result$residuals
  class(object) <- "tvsure"
  return(object)
}
##' @rdname update.tvReg
##' @return An object of class \code{tvlm}.
##' @method update tvlm
##' @export
update.tvlm <- function(object, y = NULL, ...)
{
  if(!is.null(y))
    object$y <- y
  result <- tvOLS(x = object$x, y = object$y, z = object$z,
                  bw = object$bw, est = object$est, tkernel = object$tkernel)
  object$tvcoef <- result$tvcoef
  object$fitted <- result$fitted
  object$residuals <- result$residuals
  class(object) <- "tvlm"
  return(object)
}

##' @rdname update.tvReg
##' @return An object of class \code{tvvar}.
##' @method update tvvar
##' @export
update.tvvar <- function(object, y = NULL, ...)
{
  neq <- object$neq
  obs <- object$obs
  p <- object$p
  ynames <- colnames(object$datamat)[1:neq]
  yend <- object$datamat[, 1:neq]
  rhs <- object$datamat[, -c(1:neq)]
  if(!is.null(y))
  {
    rhs <- stats::embed(y, dimension = p + 1)[, -(1:neq)]
    if(object$type == "const")
      rhs <- cbind(rhs, rep(1, obs))
    yend <- y[-c(1:p), ]
  }
  bw <- object$bw
  est <- object$est
  tkernel <- object$tkernel
  singular.ok <- object$singular.ok
  residuals <- object$resid
  fitted <- object$fitted
  equation <- object$tvcoef
  eqnames <- names(equation)
  for (i in 1:neq)
  {
    y <- yend[, i]
    results <- tvOLS(x = rhs, y = y, bw = bw[i], est = est, tkernel = tkernel,
                     singular.ok = singular.ok)
    equation[[eqnames[i]]] <- results$tvcoef
    residuals[, i] <- results$resid
    fitted[, i] <- results$fitted
  }
  object$tvcoef <- equation
  object$datamat[, 1:neq] <- yend
  object$fitted <- fitted
  object$resid <- residuals
  class(object) <- "tvvar"
  return(object)
}
