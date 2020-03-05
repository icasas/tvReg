#' Print results of functions in tvReg
#' 
#' Print some results for objects with class attribute \code{tvlm}, \code{tvar}, \code{tvvar},
#' \code{tvirf}, \code{tvsure} and \code{tvplm}.
#' 
#' These functions print a few results from the time-varying estimated coefficients
#' @param object An object used to select a method.
#' @param digits Integer, indicating the minimal number of significant digits.
#' @param ... Other parameters passed to specific methods.
#' @seealso \code{\link{plot.tvlm}}, \code{\link{plot.tvvar}}, \code{\link{plot.tvvar}},
#' \code{\link{plot.tvirf}},\code{\link{plot.tvsure}}
#' @aliases summary summary.tvlm
#' @rdname summary.tvReg
#' @method summary tvlm
#' @export

summary.tvlm <- function(object,  digits = max(3, getOption("digits") - 3), ... ) 
{
   cat("\nCall: \n")
   print(object$call)
   cat("\nClass: ", class(object), "\n")
   result <- object$coefficients
   bw <- round(object$bw, digits = digits)
   lower <- object$Lower
   upper <- object$Upper
   level <- object$level * 100
   obs <- object$obs
   text1 <- "\nSummary of time-varying estimated coefficients:"
   cat(text1, "\n")
   row <- paste(rep("=", nchar(text1)), collapse = "")
   cat(row, "\n")
   print(apply(result, 2, summary), digits = digits)
   if(!is.null(lower))
   {
     text1 <- paste("\nLOWER (", level, "%) confidence interval:", sep ="")
     cat(text1, "\n")
     print(apply(lower, 2, summary), digits = digits)
     text1 <- paste("\nUPPER (", level, "%) confidence interval:", sep ="")
     cat(text1, "\n")
     print(apply(upper, 2, summary), digits = digits)
   }
   cat("\nBandwidth: ", bw)
   SSR <- sum(object$residuals^2)
   SST <- sum((object$y - mean(object$y))^2)
   R2 <- round(1 - SSR/SST, digits = digits)
   cat("\nPseudo R-squared: ", R2, "\n\n")
   invisible(object)
}

#' @rdname summary.tvReg
#' @method summary tvar
#' @export 
summary.tvar <- summary.tvlm

#' @rdname summary.tvReg
#' @method summary tvplm
#' @export 
summary.tvplm <- summary.tvlm

#' @aliases summary summary.tvsure
#' @rdname summary.tvReg
#' @method summary tvsure
#' @export
summary.tvsure <- function (object,  digits = max(3, getOption("digits") - 3), ...)
{
  result <- object$coefficients
  bw <- round(object$bw, digits = digits)
  if (length(bw) == 1)
    bw <- rep(bw, object$neq)
  lower <- object$Lower
  upper <- object$Upper
  level <- object$level * 100
  neq <- object$neq
  nvar <- object$nvar
  names <- names(object$x)
  for (i in 1:neq)
  {
    text1 <- paste("\nSummary of tvSURE for equation \"", 
                     names[i], "\"", sep ="")
    cat(text1, "\n")
    row <- paste(rep("=", nchar(text1)), collapse = "")
    cat(row, "\n")
    index <- (sum(nvar[-c(i:neq)]) + 1):(sum(nvar[-c(i:neq)]) + nvar[i])
    print(apply(result[, index], 2, summary), digits = digits)
    if(!is.null(lower))
    {
      text1 <- paste("\nLOWER (", level, "%) confidence interval:", sep ="")
      cat(text1, "\n")
      print(apply(lower[, index], 2, summary), digits = digits)
      text1 <- paste("\nUPPER (", level, "%) confidence interval:", sep ="")
      cat(text1, "\n")
      print(apply(upper[, index], 2, summary), digits = digits)
    }
    cat("\nBandwidth: ", bw[i])
    SSR <- sum(object$residuals[, i]^2)
    SST <- sum((object$y[, i]- mean(object$y[, i]))^2)
    R2 <- round(1 - SSR/SST, digits = digits)
    cat("\nPseudo R-squared: ", R2, "\n\n")
  }
  invisible(object)
}

#' @aliases summary summary.tvvar
#' @rdname summary.tvReg
#' @method summary tvvar
#' @export
summary.tvvar <- function (object,  digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall: \n")
  print(object$call)
  result <- object$coefficients
  bw <- round(object$bw, digits = digits)
  if (length(bw) == 1)
    bw <- rep(bw, object$neq)
  neq <- object$neq
  nvar <- NCOL(object$x) - neq
  names <- colnames(object$y.orig)
  for (i in 1:neq)
  {
    text1 <- paste("\nSummary of tvVAR for equation ", 
                   names[i], ":", sep = "")
    cat(text1, "\n")
    row <- paste(rep("=", nchar(text1)), collapse = "")
    cat(row, "\n")
    print(apply(result[[i]], 2, summary), digits = digits)
    cat("\nBandwidth: ", bw[i])
    SSR <- sum(object$residuals[, i]^2)
    SST <- sum((object$y[, i]- mean(object$y[, i]))^2)
    R2 <- round(1 - SSR/SST, digits = digits)
    cat("\nPseudo R-squared: ", R2, "\n\n")
  }
  invisible(object)
}

#' @aliases summary summary.tvirf
#' @rdname summary.tvReg
#' @method summary tvirf
#' @export
summary.tvirf <-function(object, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall: \n")
    print(object$call)
    irf <- object$irf
    impulse <- object$impulse
    response <- object$response
    lower <- object$lower
    upper <- object$upper
    text1 <- NULL
    if (object$cumulative)
      text1 <- paste (text1, "Cumulative ", sep ="")
    if (object$ortho)
      text1 <- text1 <- paste(text1, "Orthogonal ", sep ="")
    text1 <- paste (text1, "tvIRF Estimation Results:", sep ="")
    cat(paste("\n", text1, "\n", sep = ""))
    row <- paste(rep("=", nchar(text1)), collapse = "")
    cat(row, "\n")
    cat("\n")
    for (i in 1:length(object$impulse))
    {
      for (j in 1:length(object$response))
      {
        text1 <- paste("\nSummary of tvIRF for impulse \"", 
                       impulse[i], "\" and response \"", 
                       response[j], "\"", sep ="")
        cat(text1, "\n")
        row <- paste(rep("=", nchar(text1)), collapse = "")
        cat(row, "\n")
        row <- apply(irf[[i]][, j, ], 2, summary)
        colnames(row) <- paste("horizon_", 0:object$n.ahead, sep ="")
        print(row, digits = digits)
      }
    }
    cat("\nBandwidth for Sigma estimation: ", 
        round(object$bw.cov, digits = digits), "\n")
    invisible(object)
}


