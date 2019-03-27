#' Print results of functions in tvReg
#' 
#' Print some results for objects with class attribute \code{tvlm}, \code{tvar}, \code{tvvar},
#' \code{tvirf}, \code{tvsure}.
#' 
#' These functions print a few results from the time-varying estimated coefficients
#' @param x An x used to select a method.
#' @param digits An integer, indicating the minimal number of significant digits.
#' @param ... Other parameters passed to specific methods.
#' @seealso \code{\link{plot.tvlm}}, \code{\link{plot.tvvar}}, \code{\link{plot.tvvar}},
#' \code{\link{plot.tvirf}},\code{\link{plot.tvsure}}
#' @rdname print.tvReg
#' @method print tvlm
#' @export

print.tvlm <- function(x,  digits = max(3, getOption("digits") - 3), ... ) 
{
   cat("\nClass: ", class(x),"\n")
   result <- x$tvcoef
   bw <- round(x$bw, digits = digits)
   lower <- x$Lower
   upper <- x$Upper
   level <- x$level * 100
   text1 <- "\nMean of coefficient estimates:"
   cat(text1, "\n")
   row <- paste(rep("=", nchar(text1)), collapse = "")
   cat(row, "\n")
   print(apply(result, 2, mean), digits = digits)
   if(!is.null(lower))
   {
    cat(paste("\nLOWER (",level, "%):", sep =""), "\n")
    print(apply(lower, 2, mean), digits = digits)
    cat(paste("\nUPPER (",level,"%):", sep =""), "\n")
    print(apply(upper, 2, mean), digits = digits)
   }
   cat("\nBandwidth: ", bw, "\n\n")
   invisible(x)
}
#' @rdname print.tvReg
#' @method print tvar
#' @export 
print.tvar <- print.tvlm

#' @inheritParams print.tvlm
#' @rdname print.tvReg
#' @method print tvsure
#' @export
print.tvsure <- function (x,  digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nClass: ", class(x),"\n")
  result <- x$tvcoef
  neq <- x$neq
  nvar <- x$nvar
  names <- names(x$x)
  method <- x$method
  bw <- round(x$bw, digits = digits)
  if (length(bw) == 1)
    bw <- rep (bw, x$neq)
  for (i in 1:neq)
  {
    text1 <- paste("\nMean of TV-SURE coefficient estimates for equation \"", 
                     names[i], "\":", sep ="")
      cat(text1, "\n")
      row <- paste(rep("=", nchar(text1)), collapse = "")
      cat(row, "\n")
      print(apply(result[, (sum(nvar[-c(i:neq)]) + 1):(sum(nvar[-c(i:neq)]) + nvar[i])], 
                  2, mean), digits = digits)
      cat("\nBandwidth: ", bw[i])
      if (method == "tvFGLS")
        cat("\t\t Covariance bandwidth: ", x$bw.cov)
      cat("\n\n")
  }
  invisible(x)
}

#' @inheritParams print.tvlm
#' @rdname print.tvReg
#' @method print tvvar
#' @export
print.tvvar <- function (x,  digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nClass: ", class(x),"\n")
  result <- x$tvcoef
  neq <- x$neq
  nvar <- ncol(x$datamat) - neq
  names <- colnames(x$y)
  bw <- round(x$bw, digits = digits)
  if (length(bw) == 1)
    bw <- rep (bw, x$neq)
  text1 <- "tvVAR Estimation Results:"
  cat(paste("\n", text1, "\n", sep = ""))
  row <- paste(rep("=", nchar(text1)), collapse = "")
  cat(row, "\n")
  cat("\n")
  for (i in 1:neq)
  {
    text1 <- paste("Mean of TV-VAR coefficient estimates for equation \"", 
                   names[i], "\":", sep = "")
    cat(text1, "\n")
    row <- paste(rep("=", nchar(text1)), collapse = "")
    cat(row, "\n")
    print(apply(result[[i]], 2, mean), digits = digits)
    cat("\nBandwidth: ", bw[i], "\n\n")
    cat("\n\n")
  }
  invisible(x)
}

#' @inheritParams print.tvlm
#' @rdname print.tvReg
#' @method print tvirf
#' @export
print.tvirf <-function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall: \n")
    print(x$call)
    irf <- x$irf
    impulse <- x$impulse
    response <- x$response
    lower <- x$lower
    upper <- x$upper
    text1 <- NULL
    if (x$cumulative)
      text1 <- paste (text1, "Cumulative ", sep ="")
    if (x$ortho)
      text1 <- text1 <- paste(text1, "Orthogonal ", sep ="")
    text1 <- paste (text1, "tvIRF Estimation Results:", sep ="")
    cat(paste("\n", text1, "\n", sep = ""))
    row <- paste(rep("=", nchar(text1)), collapse = "")
    cat(row, "\n")
    cat("\n")
    for (i in 1:length(x$impulse))
    {
      for (j in 1:length(x$response))
      {
        text1 <- paste("\nSummary of tvIRF for impulse \"", 
                       impulse[i], "\" and response \"", 
                       response[j], "\"", sep ="")
        cat(text1, "\n")
        row <- paste(rep("=", nchar(text1)), collapse = "")
        cat(row, "\n")
        row <- apply(irf[[i]][, j, ], 2, summary)
        colnames(row) <- paste("horizon_", 0:x$n.ahead, sep ="")
        print(row, digits = digits)
      }
    }
    cat("\nBandwidth for Sigma estimation: ", 
        round(x$bw.cov, digits = digits), "\n")
    invisible(x)
}


