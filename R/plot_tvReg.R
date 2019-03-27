#' Plot Methods for Objects in tvReg
#'
#' Plot methods for objects with class attribute \code{tvlm}, \code{tvar}, \code{tvvar},
#' \code{tvirf}, \code{tvsure}.
#' @rdname plot.tvReg
#' @method plot tvsure
#' @param x An x used to select a method.
#' @param ... Other parameters passed to specific methods.
#' @param eqs A vector of integers. Equation(s) number(s) of the coefficients to be plotted.
#' @param vars A vector of integers. Variable number(s) of the coefficients to be plotted.
#' @param plot.type	Character, if multiple all plots are drawn in a single device,
#' otherwise the plots are shown consecutively.
#' @seealso \code{\link{tvLM}}, \code{\link{tvAR}}, \code{\link{tvVAR}},
#' \code{\link{tvSURE}}
#'
#' @export
#'
plot.tvsure <- function(x, eqs = NULL, vars = NULL, 
                        plot.type = c("multiple", "single") , ...)
{
  if (class(x) != "tvsure")
    stop("\nPlot not implemented for this class.\n")
  if(!is.null(vars) & any(vars <= 0))
    stop("\nInvalid number in 'vars'\n")
  if(!is.null(eqs) & any(eqs <= 0))
    stop("\nInvalid number in 'eqs'\n")
  if(!is.null(eqs) & any(eqs > x$neq))
    stop("\nInvalid number in 'eqs'\n")
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  if (!any(plot.type %in% c("multiple", "single")))
    stop("\nParameter plot.type only takes values \"multiple\" and \"single\". \n")
  tvcoef <- x$tvcoef
  if(is.null(tvcoef))
    stop("\nThe time-varying coefficients array is NULL. \n")
  Lower <- x$Lower
  Upper <- x$Upper
  nvar <- x$nvar
  neq <- x$neq
  obs <- x$obs
  z <- x$z
  eq.names <- names(x$x)
  if (is.null(eqs))
    eqs <- 1:neq
  if (max(eqs) > neq)
    stop("\nSome of the equation chosen are not part of the model.\n")
  x.axis <- 1:obs
  sub <- "t"
  if(!is.null(z))
  {
    sort.index <- sort.int(z, index.return = TRUE)$ix
    x.axis <- z[sort.index]
    tvcoef <- tvcoef[sort.index, , , drop = FALSE]
    sub <- expression(z[t])
    if(!is.null(lower))
    {
      Lower <- Lower[sort.index, , drop = FALSE]
      Upper <- Upper[sort.index, , drop = FALSE]
    }
  }
  if (x$level != 0)
    sub <- paste0(sub, "\n", (x$level) * 100, "% Bootstrap CI, ",
                 x$runs, " runs")
  for (i in eqs)
  {
    var.names <- colnames(x$x[[i]])
    if(any(vars > nvar[i]))
      warning("\nSome variables you want to plot do not exist in equation ", i, "\n")
    if(!is.null(vars))
      var.names <- var.names[vars]
    plotvars <- which(colnames(x$x[[i]]) %in% var.names)
    coef <- tvcoef[, (sum(nvar[1:i])-nvar[i]+1):sum(nvar[1:i]), drop = FALSE]
    lower <- Lower[, (sum(nvar[1:i])-nvar[i]+1):sum(nvar[1:i]), drop = FALSE]
    upper <- Upper[, (sum(nvar[1:i])-nvar[i]+1):sum(nvar[1:i]), drop = FALSE]
    coef <- coef[, plotvars, drop = FALSE]
    lower <- lower[, plotvars, drop = FALSE]
    upper <- upper[, plotvars, drop = FALSE]
    if(length(plotvars) == 1)
      plot.type <- "single"
    if(any(plot.type == "multiple"))
    {
      nplots <- ceiling (length(plotvars)/3)
      graphics::par(mfrow = c(min(3, length(plotvars)), 1), 
                    mar = c(0, 4, 0, 1), oma = c(6, 4, 3, 1))
      count <- 1
      while (count <= nplots)
      {
        for (j in ((count-1)*3 + 1):(min(count*3, length(plotvars))))
        {
          ylim <- range(coef[, j], upper[, j], lower[, j])
          graphics::plot(x.axis, coef[, j], axes = FALSE, type = "l", ylab = var.names[j],
                         ylim = ylim, ...)
          graphics::axis(2, at = pretty(ylim)[-1])
          if(!is.null(lower))
          {
            graphics::polygon(c(rev(x.axis), x.axis), c(rev(upper[, j]), lower[, j]),
                              col = "grey80", border = NA, fillOddEven = TRUE)
            graphics::lines(x.axis, coef[, j])
          }
          graphics::abline(h = 0, col = 2, lty = 2)
          graphics::box()
        }
        graphics::axis(1, at = pretty(x.axis)[-1])
        graphics::mtext(eq.names[i], 3, line = 1, outer = TRUE, ...)
        graphics::mtext(sub, 1, line = 4, outer = FALSE, ...)
        count <- count + 1
      }
      graphics::par(ask = TRUE)
    }
    else
    {
      for (j in 1:length(plotvars))
      {
        graphics::par(mfrow = c(1, 1), mar = c(5, 5, 0, 1), oma = c(2, 0, 6, 0))
        ylim <- range(coef[, j], upper[, j], lower[, j])
        graphics::plot(x.axis, coef[, j], axes = FALSE, type = "l", 
                       ylab = var.names[j], ylim = ylim, xlab ="", ...)
        graphics::axis(2, at = pretty(ylim)[-1])
        if(!is.null(lower))
        {
          graphics::polygon(c(rev(x.axis), x.axis), c(rev(upper[, j]), lower[, j]),
                            col = "grey80", border = NA, fillOddEven = TRUE)
          graphics::lines(x.axis, coef[, j])
        }
        graphics::abline(h = 0, col = 2, lty = 2)
        graphics::box()
        graphics::axis(1, at = pretty(x.axis)[-1])
        graphics::mtext(eq.names[i], 3, line = 2, outer = FALSE, ...)
        graphics::mtext(sub, 1, line = 3, outer = FALSE, ...)
        graphics::par(ask = TRUE)
      }
    }
  }
}
#' @rdname plot.tvReg
#' @method plot tvlm
#' @export
#'
plot.tvlm <- function(x, ...)
{
  if (!any(class(x) %in% c("tvlm", "tvar")))
    stop("\nPlot not implemented for this class.\n")
  .univariatePlot (x)
}

#' @rdname plot.tvReg
#' @method plot tvar
#' @export 
plot.tvar <- plot.tvlm

#' @name tvReg-internals
#' @aliases .univariatePlot
#' @title tvReg internal and secondary functions
#' @keywords internal
.univariatePlot <-function(x, ...)
{
  tvcoef <- x$tvcoef
  if(is.null(tvcoef))
    stop("\nThe time-varying coefficients matrix is NULL. \n")
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  lower <- x$Lower
  upper <- x$Upper
  nvar <- NCOL(tvcoef)
  obs <- NROW(tvcoef)
  var.names <- colnames(tvcoef)
  z <- x$z
  sub <- "t"
  x.axis <- 1:obs
  if(!is.null(z))
  {
    sort.index <- sort.int(z, index.return = TRUE)$ix
    x.axis <- z[sort.index]
    tvcoef <- x$tvcoef[sort.index, , drop = FALSE]
    xlabel <- "z"
    if(!is.null(lower))
    {
      lower <- lower[sort.index, , drop = FALSE]
      upper <- upper[sort.index, , drop = FALSE]
    }
    sub <- expression(z[t])
  }
  if (x$level != 0)
    sub <- paste0(sub, "\n", (x$level) * 100, "% Bootstrap CI, ",
                 x$runs, " runs")
  graphics::par(mfrow = c(1, 1), 
                mar = c(4, 4, 2, 1), oma = c(0, 0, 0, 0))
  for ( j in 1:nvar)
  {
    ylim <- range(tvcoef[, j])
    if(!is.null (lower))
      ylim <- range(ylim, lower[, j], upper[, j])
    graphics::plot(x.axis, tvcoef[, j], xlab = "", ylab = var.names[j], 
                   type = "l", ylim = ylim, axes = FALSE, ...)
    if(!is.null(lower))
    {
      graphics::polygon(c(rev(x.axis), x.axis), c(rev(upper[, j]), lower[, j]),
                        col = "grey80", border = NA, fillOddEven = TRUE)
      graphics::lines(x.axis, tvcoef[, j])
      
    }
    graphics::axis(2, at = pretty(ylim)[-1])
    graphics::axis(1, at = pretty(x.axis)[-1])
    graphics::box()
    graphics::mtext(sub, 1, line = 3, outer = FALSE, ...)
    if(nvar > 1)
      graphics::par(ask = TRUE)
  }
}

#' @rdname plot.tvReg
#' @method plot tvvar
#' @export
#'
plot.tvvar <- function(x, ...)
{
  if (class(x) != "tvvar")
    stop("\nPlot not implemented for this class.\n")
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  fitted <- x$fitted
  residuals <- x$resid
  neq <- x$neq
  obs <- x$obs
  y <- x$y
  var.names <- colnames(y)
  graphics::par(mar = c(2, 3, 2, 1), mfrow = c(2,1), cex = 1)
  for ( j in 1:neq)
  {
    ylim <- range(y[, j], fitted[, j])
    graphics::plot(1:obs, y[, j], ylim = ylim, xlab = "", ylab  = "",
                   main = paste("Diagram of fit for ", var.names[j], sep = ""),
                   pch = 20, cex = 0.5, yaxt = "n")
    graphics::lines(1:obs, fitted[, j], col = 2)
    graphics::axis(2, at = pretty(ylim)[-1])
    ylim <- range(residuals[, j])
    graphics::plot(1:obs, residuals[, j], xlab = "", ylab = "", type = "l",
                   main = paste("Diagram of residuals for ", var.names[j], 
                                sep = ""), yaxt ="n")
    graphics::axis(2, at = pretty(ylim)[-1])
    graphics::abline(h = 0, lty = 2)
    if(neq > 1)
      graphics::par(ask = TRUE)
  }
}

#' @rdname plot.tvReg
#' @method plot tvirf
#' @param impulse	Character  vector (optional) of the impulses, default is all variables.
#' @param response Character vector (optional) of the responses, default is all variables.
#' @param obs.index  Scalar (optional), the time at which the impulse response is plotted.
#' If left NULL, the mean over the whole period is plotted (this values should be similar to
#' the estimation using a non time-varying VAR method).
#' @inheritParams plot.tvsure
#' @export
#'
plot.tvirf <- function (x, obs.index = NULL, impulse = NULL, response = NULL,
                        plot.type = c("multiple", "single"), ...)
{
  if(is.null(obs.index))
  {
    cat("\nThe plot represents the mean of tvIRF over every time period. 
          Enter the row number in parameter \"obs.index\" to plot the tvIRF 
          of a particular point in time.\n")
  }
  else if (!is.null(obs.index) & length (obs.index) > 1)
    stop("\nPlease enter only one value in 'obs.index' or nothing to get the average of
         the whole period.\n")
  else if(!is.null(obs.index) & (obs.index > x$x$obs | obs.index <= 0))
    stop("\nWrong index: too large or too small")
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  plot.type <- match.arg(plot.type)
  if(is.null(impulse))
    impulse <- x$impulse
  if(is.null(response))
    response <- x$response
  if(sum(response %in% x$response)< length (response))
    stop ("\nOne or several response variables are not part of the model.\n")
  if(sum(impulse %in% x$impulse)< length (impulse))
    stop ("\nOne or several impulse variables are not part of the model.\n")
  inames <- impulse
  rnames <- response
  index <- response %in% x$response
  if(length(rnames) == 1)
    plot.type <- "single"
  if (is.null(obs.index))
    obs.index2 <- 1:x$x$obs
  else
    obs.index2 <- obs.index
  x.axis <- 1:tail(dim(x$irf[[1]]), 1)
  for (i in 1:length(impulse))
  {
    iname <- inames[i]
    main <- paste("Impulse variable: ", iname, sep="")
    if (x$cumulative)
      main <- paste(main, "(cumulative)", sep = " ")
    sub <- "horizon"
    if (x$level != 0)
      sub <- paste0(sub, "\n",(x$level) * 100, "% Bootstrap CI, ",
                     x$runs, " runs")
    ylabel = rnames
    impulses <- x$irf[[iname]][obs.index2, rnames, ,drop = FALSE]
    impulses <- apply (impulses, 2:3, mean)
    upper <- NULL
    lower <- NULL
    if (x$level != 0) 
    {
      upper <- x$Upper[[iname]][obs.index2, rnames, , drop = FALSE]
      lower <- x$Lower[[iname]][obs.index2, rnames, , drop = FALSE]
      upper <- apply (upper, 2:3, mean)
      lower <- apply (lower, 2:3, mean)
    }
    nresponse <- length(rnames)
    nplots <- ceiling (nresponse/4)
    if(nresponse == 1)
      plot.type <- "single"
    if(plot.type == "multiple")
    {
      graphics::par(mfrow = c(min(4, nresponse), 1), mar = c(0, 4, 0, 4), oma = c(6, 4, 6, 4))
      count <- 1
      while (count <= nplots)
      {
        for (j in ((count-1)*4 +1):(min(count*4, nresponse)))
        {
          rname <- rnames[j]
          ylim <- range(impulses[rname, ], upper[rname, ], lower[rname, ])
          graphics::plot(x.axis, impulses[rname, ], axes = FALSE, type = "l", ylab = ylabel[j],
                         ylim = ylim, ...)
          graphics::axis(2, at = pretty(ylim)[-1])
          if(!is.null(lower))
          {
            graphics::polygon(c(rev(x.axis), x.axis), c(rev(upper[rname, ]), lower[rname, ]),
                              col = "grey80", border = NA, fillOddEven = TRUE)
            graphics::lines(x.axis, impulses[rname, ])
          }
          graphics::abline(h = 0, col = 2, lty = 2)
          graphics::box()
        }
        count <- count + 1
        graphics::axis(1, at = x.axis, labels = c(x.axis - 1))
        graphics::mtext(main, 3, line = 2, outer = TRUE, ...)
        graphics::mtext(sub, 1, line = 4, outer = FALSE, ...)
        graphics::par(ask = TRUE)
      }
     
    }
    else if(plot.type == "single")
    {
      graphics::par(mfrow = c(1, 1), mar = c(4, 4, 4, 1), oma = c(3, 0, 2, 1))
      ylabel <- rnames
    
      for (j in 1:length(response))
      {
        rname <- rnames[j]
        ylim <- range(c(impulses[rname, ], lower[rname, ], upper[rname, ]))
        graphics::plot(x.axis, impulses[rname, ], type = "l", ylim = ylim, axes = FALSE,
                     ylab = ylabel[j], main = main, xlab = "", ...)
        if(!is.null(lower))
        {
          graphics::polygon(c(rev(x.axis), x.axis), c(rev(upper[rname, ]), lower[rname,]),
                            col = "grey80", border = NA, fillOddEven = TRUE)
          graphics::lines(x.axis, impulses[rname, ])
        }
        graphics::abline(h = 0, col = 2, lty = 2)
        graphics::axis(1, at = x.axis, labels = c(x.axis-1))
        graphics::axis(2, at = pretty(ylim)[-1])
        graphics::mtext(main, 3, line = 2, outer = TRUE, ...)
        graphics::mtext(sub, 1, line = 3, outer = FALSE, ...)
        graphics::box()
        graphics::par(ask = TRUE)
      }
    }
  }
}
