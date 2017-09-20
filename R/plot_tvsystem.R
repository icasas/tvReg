##' Plot Methods for objects in tvReg
##'
##' Plot methods for objects with class attribute \code{tvlm}, \code{tvar}, \code{tvvar},
##' \code{tvirf}, \code{tvsure}.
##' @rdname plot.tvsure
##' @method plot tvsure
##' @param x An x used to select a method.
##' @param ... Other parameters passed to specific methods.
##' @param eqs Character vector (optional) with the equation(s) number(s) or
##' equation name(s) of the coefficients to be plotted.
##' @param vars Character vector (optional) with the variable number(s) or
##' variable name(s) of the coefficients to be plotted.
##' @seealso \code{\link{tvLM}}, \code{\link{tvAR}}, \code{\link{tvVAR}},
##' \code{\link{tvSURE}}
##' @examples
##'
##' data( "Kmenta", package = "systemfit" )
##' eqDemand <- consump ~ price + income
##' eqSupply <- consump ~ price + farmPrice + trend
##' system <- list( demand = eqDemand, supply = eqSupply )
##'
##' ## tvOLS estimation
##' tvols.fit <- tvSURE(system, data = Kmenta)
##' ## 95% confidence interval using Mammen's wild bootstrap.
##' tvols.fit <- CI(tvols.fit, level = 0.95)
##' plot(tvols.fit)
##'
##' ## Time-varying Five Factor Model
##' data(FF5F)
##' ## SMALL/LoBM porfolios time-varying three factor model
##' eqNA<-NA.SMALL.LoBM ~ NA.Mkt.RF + NA.SMB+ NA.HML
##' eqJP<-JP.SMALL.LoBM ~ JP.Mkt.RF + JP.SMB+ JP.HML
##' eqAP<-AP.SMALL.LoBM ~ AP.Mkt.RF + AP.SMB+ AP.HML
##' eqEU<-EU.SMALL.LoBM ~ EU.Mkt.RF + EU.SMB+ EU.HML
##' system<-list(NorthA = eqNA, JP = eqJP, AP = eqAP, EU = eqEU)
##'
##' ## Fit a time-varying coefficients SURE model
##' ff5f.tv<-tvSURE(system, data = FF5F, method = "tvFGLS", est = "ll",
##' bw = c(0.11, 0.43, 0.49, 0.27))
##' ## 95% confidence interval using Mammen's wild bootstrap.
##' ff5f.tv <- CI(ff5f.tv, level = 0.95, runs = 30)
##'
##' ## Plot the intercepts
##' plot(ff5f.tv, vars = 1)
##'
##' @export
##'
plot.tvsure <- function(x, eqs = NULL, vars = NULL, ...)
{
  if (class(x) != "tvsure")
    stop("\nPlot not implemented for this class.\n")
  tvcoef <- x$tvcoef
  if(is.null(tvcoef))
    stop("\nThe time-varying coefficients array is NULL. \n")
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  Lower <- x$Lower
  Upper <- x$Upper
  nvar <- x$nvar
  neq <- length(x$x)
  obs <- x$obs
  z <- x$z
  eq.names <- names(x$x)
  y.names <- colnames(x$y)
  if (is.null(eqs))
    eqs <- 1:neq
  x.axis <- 1:obs
  if(!is.null(z))
  {
    sort.index <- sort.int(z, index.return = TRUE)$ix
    x.axis <- z[sort.index]
    tvcoef <- tvcoef[sort.index, , , drop = FALSE]
    x.lab <- "z"
    if(!is.null(Lower))
    {
      Lower <- Lower[sort.index, , drop = FALSE]
      Upper <- Upper[sort.index, , drop = FALSE]
    }
  }
  graphics::par(mar = c(4, 4, 4, 1))
  for (i in eqs)
  {
    var.names<-colnames(x$x[[i]])
    coef <- tvcoef[, (sum(nvar[1:i])-nvar[i]+1):sum(nvar[1:i])]
    lower <- Lower[, (sum(nvar[1:i])-nvar[i]+1):sum(nvar[1:i])]
    upper <- Upper[, (sum(nvar[1:i])-nvar[i]+1):sum(nvar[1:i])]
    if(is.null(vars))
      vars <- 1:nvar[i]
    for ( j in vars)
    {
      graphics::plot(x.axis, coef[, j], main = paste(eq.names[i], ": effect over time on ",
                                             y.names[i], sep = ""),
                     xlab = "time", ylab = var.names[j], type = "l",
                     ylim = range(coef[, j], lower[, j], upper[, j]))
      if(!is.null(lower))
      {
        graphics::polygon(c(rev(x.axis), x.axis), c(rev(upper[, j]), lower[, j]),
                          col = "grey80", border = NA, fillOddEven = TRUE)
        graphics::lines(x.axis, coef[, j])
      }
      graphics::par(ask = TRUE)
    }
  }
}

##' @rdname plot.tvsure
##' @method plot tvlm
##' @export
##'
plot.tvlm <- function(x, ...)
{
  if (!any(class(x) == "tvlm"))
    stop("\nPlot not implemented for this class.\n")
  tvcoef <- x$tvcoef
  if(is.null(tvcoef))
    stop("\nThe time-varying coefficients matrix is NULL. \n")
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  Lower <- x$Lower
  Upper <- x$Upper
  nvar <- ncol(tvcoef)
  obs <- nrow(tvcoef)
  var.names <- colnames(tvcoef)
  graphics::par(mar = c(4, 4, 2, 1))
  z <- x$z
  x.lab <- "t"
  x.axis <- 1:obs
  if(!is.null(z))
  {
    sort.index <- sort.int(z, index.return = TRUE)$ix
    x.axis <- z[sort.index]
    tvcoef <- x$tvcoef[sort.index, , drop = FALSE]
    x.lab <- "z"
    if(!is.null(Lower))
    {
      Lower <- Lower[sort.index, , drop = FALSE]
      Upper <- Upper[sort.index, , drop = FALSE]
    }
  }
  for ( j in 1:nvar)
  {
    ylim <- range(tvcoef[, j])
    if(!is.null (Lower))
      ylim <- range(ylim, Lower[, j], Upper[, j])
    graphics::plot(x.axis, tvcoef[, j], xlab = x.lab, ylab = var.names[j], type = "l",
                   ylim = ylim)
    if(!is.null(Lower))
    {
      graphics::polygon(c(rev(x.axis), x.axis), c(rev(Upper[, j]), Lower[, j]),
                        col = "grey80", border = NA, fillOddEven = TRUE)
      graphics::lines(x.axis, tvcoef[, j])
    }
    if(nvar > 1)
      graphics::par(ask = TRUE)
  }
}

##' @rdname plot.tvsure
##' @method plot tvvar
##' @export
##'
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
  y <- x$datamat[, 1:neq]
  var.names <- colnames(y)
  graphics::par(mar = c(4, 4, 2, 1), mfrow = c(2,1))
  for ( j in 1:neq)
  {
    ylim <- range(y[, j], fitted[, j])
    graphics::plot(1:obs, y[, j], ylim = ylim, xlab = "", ylab  = "",
                   main = paste("Diagram of fit for ", var.names[j], sep = ""),
                   pch = 20, cex = 0.5)
    graphics::lines(1:obs, fitted[, j], col = 2)
    graphics::plot(1:obs, residuals[, j], xlab = "", ylab = "", type = "l",
                   main = paste("Diagram of residuals for ", var.names[j], sep = ""))
    graphics::abline(h = 0, lty = 2)
    if(neq > 1)
      graphics::par(ask = TRUE)
  }
}


##' @rdname plot.tvsure
##' @method plot tvirf
##' @param impulse	Character  vector (optional) of the impulses, default is all variables.
##' @param response Character vector (optional) of the responses, default is all variables.
##' @param adj.mtext	Adjustment for mtext(), only applicable if plot.type = "multiple".
##' @param obs.index  Scalar (optional), the time at which the impulse response is plotted.
##' If left NULL, the mean over the whole period is plotted (this values should be similar to
##' the estimation using a non time-varying VAR method).
##' @param main Character vector, the titles of the plot.
##' @param mar.multi	Setting of margins, if plot.type = "multiple".
##' @param names	Character vector (optional), the variables names to be plotted.
##' If left NULL, all variables are plotted.
##' @param sub Character, sub title in plot.
##' @param nc Integer, number of columns for multiple plot.
##' @param oma.multi	Setting of margins, if plot.type = "multiple".
##' @param padj.mtext Adjustment for mtext(), only applicable if plot.type = "multiple".
##' @param plot.type	Character, if multiple all plots are drawn in a single device,
##' otherwise the plots are shown consecutively.
##' @param xlab	Character vector signifying the labels for the x-axis.
##' @param ylab Character vector signifying the labels for the y-axis.
##'
##' @examples
##' ## Inflation rate, unemployment rate and treasury bill interest rate for the US,
##' ## as used by Primiceri (2005).
##'
##' data(usmacro, package = "bvarsv")
##'
##' ##Estimate a time-varying coefficients vector autoregressive and its
##' ##impulse response function
##' model.tvVAR <- tvVAR(usmacro, p = 6, type = "const")
##' model.tvIRF <- tvIRF(model.tvVAR)
##'
##' ## Plot residuals and fitted values of tvVAR
##'
##' ##Obtain 95% confidence interval of the impulse response function
##' model.tvIRF <- CI(model.tvIRF)
##'
##' ##plot the mean tvIRF over the whole period
##' plot(model.tvIRF)
##'
##' ##plot the tvIRF at time 50
##' plot(model.tvIRF, obs.index = 50)
##'
##' ##plot the effect of a shock in inflation on unemployment
##' plot(model.tvIRF, impulse = "inf", response = "une")
##'
##' @export
##'
plot.tvirf <- function (x, obs.index = NULL, impulse = NULL, response = NULL,
                        plot.type = c("multiple", "single"),
                        names = NULL, main = NULL, sub = NULL, ylab = NULL,
                        xlab = NULL, nc, mar.multi = c(0, 4, 0, 4),
                        oma.multi = c(6, 4, 6, 4), adj.mtext = NA, padj.mtext = NA, ...)
{
  if(is.null(obs.index))
  {
    cat("\nThe mean of tvIRF over all time period will be plotted. Enter a row number in obs.index
        to plot a particular point in time\n")
  }
  else if (!is.null(obs.index) & length (obs.index) > 1)
    stop("\nPlease enter only one value in 'obs.index' or nothing to get the average of
         the whole period.\n")
  else if(!is.null(obs.index) & (obs.index > x$x$obs | obs.index <= 0))
    stop("\nWrong index: too large or too small")
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  plot.type <- match.arg(plot.type)
  inames <- impulse
  rnames <- response
  if(length(rnames) == 1)
    plot.type <- "single"
  if (is.null (impulse))
    inames <- x$impulse
  if (is.null(response))
    rnames <- x$response
  if (is.null(names)) {
    names <- inames
  }
  else {
    names <- as.character(names)
    if (!(all(names %in% inames))) {
      warning("\nInvalid variable name(s) supplied, using first variable.\n")
      inames <- inames[1]
    }
    else {
      inames <- names
    }
  }
  nvi <- length(inames)
  nvr <- length(rnames)
  dataplot <- function(x, obs.index, iname, rnames)
  {
    if (is.null(obs.index))
      obs.index2 <- 1:x$x$obs
    else
      obs.index2 <- obs.index
    impulses <- x$irf[[iname]][obs.index2, rnames, ,drop = FALSE]
    impulses <- apply (impulses, 2:3, mean)
    upper <- NULL
    lower <- NULL
    if (x$level != 0) {
      upper <- x$Upper[[iname]][obs.index2, rnames, , drop = FALSE]
      lower <- x$Lower[[iname]][obs.index2, rnames, , drop = FALSE]
      upper <- apply (upper, 2:3, mean)
      lower <- apply (lower, 2:3, mean)
    }
    text1 <- paste("Impulse variable: ", iname, sep="")
    if (x$cumulative)
      text1 <- paste(text1, "(cumulative)", sep = " ")
    text2 <- ""
    if (x$level != 0)
      text2 <- paste((x$level) * 100, "% Bootstrap CI, ",
                     x$runs, "runs")
    result <- list(impulses = impulses, upper = upper, lower = lower,
                   obs.index = obs.index, text1 = text1, text2 = text2)
    return(result)
  }
  plot.single <- function(dp, iname, rname, ...)
  {
    x <- dp$impulses
    upper <- dp$upper
    lower <- dp$lower
    ifelse(is.null(main), main <- dp$text1, main <- main)
    ifelse(is.null(sub), sub <- dp$text2, sub <- sub)
    if(length(rname) == 1)
      x.axis <- 1:length(x)
    else
      x.axis <- 1:ncol(x)
    ifelse(is.null(ylab), ylabel <- rname, ylabel <- ylab)
    ifelse(is.null(xlab), xlabel <- "", xlabel <- xlab)
    ylim <- range(c(x, lower, upper))
    graphics::plot(x.axis, x, type = "l", ylim = ylim, axes = FALSE,
         ylab = paste(ylabel), xlab = paste(xlab), ...)
    graphics::title(main = main, sub = sub, ...)
    graphics::axis(1, at = x.axis, labels = c(x.axis-1))
    graphics::axis(2, ...)
    graphics::box()
    if(!is.null(lower))
    {
      graphics::polygon(c(rev(x.axis), x.axis), c(rev(upper), lower),
                        col = "grey80", border = NA, fillOddEven = TRUE)
      graphics::lines(x.axis, x)
    }
    graphics::abline(h = 0, col = 2, lty = 2)
  }
  plot.multiple <- function(dp, nc = nc, ...)
  {
    x <- dp$impulses
    upper <- dp$upper
    lower <- dp$lower
    x.axis <- 1:ncol(x)
    nvr <- nrow(x)
    ifelse(is.null(main), main <- dp$text1, main <- main)
    ifelse(is.null(sub), sub <- dp$text2, sub <- sub)
    if (missing(nc))
      nc <- ifelse(nvr > 4, 2, 1)
    nr <- ceiling(nvr/nc)
    graphics::par(mfrow = c(nr, nc), mar = mar.multi, oma = oma.multi)
    if (nr > 1)
    {
      for (j in 1:(nvr - nc))
      {
        ifelse(is.null(ylab), ylabel <- rownames(x)[j],
               ylabel <- ylab)
        ylim <- range(x[j, ], upper[j, ], lower[j,])
        graphics::plot(x.axis, x[j, ], axes = FALSE, type = "l", ylab = ylabel,
             ylim = ylim, ...)
        graphics::axis(2, at = pretty(ylim)[-1])
        if(!is.null(lower))
        {
          graphics::polygon(c(rev(x.axis), x.axis), c(rev(upper[j, ]), lower[j, ]),
                            col = "grey80", border = NA, fillOddEven = TRUE)
          graphics::lines(x.axis, x[j, ])
        }
        graphics::abline(h = 0, col = 2, lty = 2)
        graphics::box()
      }
      for (j in (nvr - nc + 1):nvr)
      {
        ifelse(is.null(ylab), ylabel <- rownames(x)[j],
               ylabel <- ylab)
        ylim <- range(x[j, ], upper[j, ], lower[j,])
        graphics::plot(x.axis, x[j, ], axes = FALSE, type = "l", ylab = ylabel,
             ylim = ylim, ...)
        graphics::axis(2, at = pretty(ylim)[-1])
        graphics::axis(1, at = x.axis, labels = c(x.axis - 1))
        if(!is.null(lower))
        {
          graphics::polygon(c(rev(x.axis), x.axis), c(rev(upper[j, ]), lower[j, ]),
                            col = "grey80", border = NA, fillOddEven = TRUE)
          graphics::lines(x.axis, x[j, ])
        }
        graphics::box()
        graphics::abline(h = 0, col = 2, lty = 2)
      }
      graphics::mtext(main, 3, line = 2, outer = TRUE, adj = adj.mtext,
            padj = padj.mtext, ...)
      graphics::mtext(sub, 1, line = 4, outer = TRUE, adj = adj.mtext,
            padj = padj.mtext, ...)
    }
    else {
      for (j in 1:nvr)
      {
        ifelse(is.null(ylab), ylabel <- rownames(x)[j],
               ylabel <- ylab)
        graphics::plot(x.axis, x[j, ], type = "l", ylab = ylabel, ylim = ylim, ...)
        if(!is.null(lower))
        {
          graphics::polygon(c(rev(x.axis), x.axis), c(rev(upper[j, ]), lower[j, ]),
                            col = "grey80", border = NA, fillOddEven = TRUE)
          graphics::lines(x.axis, x[j, ])
        }
        graphics::abline(h = 0, col = 2, lty = 2)
      }
      graphics::mtext(main, 3, line = 2, outer = TRUE, adj = adj.mtext,
            padj = padj.mtext, ...)
      graphics::mtext(sub, 1, line = 4, outer = TRUE, adj = adj.mtext,
            padj = padj.mtext, ...)
    }
  }
  if (plot.type == "single")
  {
    for (i in 1:nvi)
    {
      dp <- dataplot(x, obs.index = obs.index, iname = inames[i], rnames = rnames)
      for (j in 1:nvr)
      {
        plot.single(dp, iname = inames[i], rname = rnames[j],
                    ...)
        if (nvr > 1)
          graphics::par(ask = TRUE)
      }
      if (nvi > 1)
        graphics::par(ask = TRUE)
    }
  }
  if (plot.type == "multiple")
  {
    for (i in 1:nvi)
    {
      dp <- dataplot(x, obs.index = obs.index, iname = inames[i], rnames = rnames)
      plot.multiple(dp, nc = nc, ...)
      if (nvi > 1)
        graphics::par(ask = TRUE)
    }
  }
}
