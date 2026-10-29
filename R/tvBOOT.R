# @rdname tvReg-internals
#' @keywords internal
#'
.tvboot<-function (x , ...) UseMethod(".tvboot", x)

# @rdname tvReg-internals
#' @method .tvboot default
#' @keywords internal
#'
.tvboot.default<-function (x, runs = 100, tboot = c("wild", "wild2", "lbwb"), ...)
{
  if (!(any(class(x) %in% c("tvlm", "tvar"))))
  {
    stop("Bootstrap confidence intervals not implemented for this class.\n")
  }
  obs <- length(x$y)
  nvar <- NCOL (x$x)
  beta.pilot <- NULL
  tboot <- match.arg(tboot)
  if(tboot == "lbwb")
  {
    h <- x$bw
    l <- ceiling(4.5 * (obs * h)^(1/4))
    model.pilot <- update(x, bw = 2*(h^(5/9)))
    resorig <- c(scale(residuals(model.pilot), scale = FALSE))
    fitted <- model.pilot$fitted
    beta.pilot <- model.pilot$coefficients
  }
  else
  {
    resorig <- c(scale(x$residuals, scale = FALSE))
    fitted <- x$fitted
    prob <- c(0.7236067977499789360962267892318777740001678466796875,
              0.2763932022500210639037732107681222259998321533203125)
    binar <- c((1-sqrt(5))*0.5, (1+sqrt(5))*0.5)
  }
  yboot <- matrix(NA, nrow = obs, ncol = runs)

  for (i in 1:runs) 
  {
    if (tboot == "wild")
      resid <- resorig*sample(binar, obs, replace = TRUE, prob = prob)
    else if (tboot == "wild2")
      resid <- resorig*stats::rnorm(obs)
    else if (tboot == "lbwb")
    {
      obs.blocks <- ceiling(obs / l)
      u.star <- numeric(obs.blocks * l)
      for (k in 1:obs.blocks)
      {
        target <- (k - 1) * l + 1
        l.index <- max(1, target - l)
        u.index <- min(obs - l + 1, target + l)
        pi.i <- sample(l.index:u.index, 1)
        
        block.selected <- resorig[pi.i:(pi.i + l - 1)]
        fill.start <- (k - 1) * l + 1
        len_to_fill <- min(l, obs - fill.start + 1)
        
        # Apply local normal multiplier over the selected residual block
        u.star[fill.start:(fill.start + len_to_fill - 1)] <- stats::rnorm(1) * block.selected[1:len_to_fill]
      }
      resid <- u.star[1:obs]
    }
    yboot[, i] <- fitted + resid
  }
  return(list(BOOT = .tvLM.ci(x, yboot), beta.pilot = beta.pilot))
    
}

# @rdname tvReg-internals
#' @method .tvboot tvplm
#' @keywords internal
#'
.tvboot.tvplm <- function (x , runs = 100, tboot = "wild", ...)
{
  obs <- x$obs
  neq <- x$neq
  resorig <- matrix(x$residuals, nrow = obs, ncol = neq)
  resorig <- apply(resorig, 2, scale, scale = FALSE)
  prob <- c(0.7236067977499789360962267892318777740001678466796875,
            0.2763932022500210639037732107681222259998321533203125)
  binar <- c((1-sqrt(5))*0.5, (1+sqrt(5))*0.5)
  fitted <- x$fitted
  yboot <- matrix(NA, nrow = obs*neq, ncol = runs)
  for (i in 1:runs)
  {
    if (tboot == "wild")
      resid <- resorig*sample(binar, obs, replace = TRUE, prob = prob)
    else if (tboot == "wild2")
      resid <- resorig*stats::rnorm(obs*neq)
    yboot[, i] <- fitted + as.numeric(resid)
  }
  return(.tvPLM.ci(x, yboot))
}

# @rdname tvReg-internals
#' @method .tvboot tvsure
#' @keywords internal
#'
.tvboot.tvsure<-function (x, runs = 100, tboot = "wild", ...)
{
  obs <- x$obs
  nvar<- x$nvar
  neq <- x$neq
  resorig <- scale(x$residuals, scale = FALSE)
  fitted <- x$fitted
  prob <- c(0.7236067977499789360962267892318777740001678466796875,
            0.2763932022500210639037732107681222259998321533203125)
  binar <- c((1-sqrt(5))*0.5, (1+sqrt(5))*0.5)
  yboot <- vector("list", runs)
  for (i in 1:runs) {
    if (tboot == "wild")
      resid <- resorig*sample(binar, obs, replace = TRUE, prob = prob)
    else if (tboot == "wild2")
      resid <- resorig*stats::rnorm(obs*neq)
    yboot[[i]] <- fitted + resid
  }
  return(.tvSURE.ci(x, yboot))
}

# @rdname tvReg-internals
#' @method .tvboot tvirf
#' @keywords internal
.tvboot.tvirf<-function (x, runs = 100, tboot = c("wild", "wild2", "lbwb"), ...)
{
  ortho <- x$ortho
  cumulative <- x$cumulative
  impulse <- x$impulse
  response <- x$response
  bw.cov <- x$bw.cov
  ortho.cov <- x$ortho.cov
  n.ahead <- x$n.ahead
  unit.shock <- x$unit.shock
  tboot <- match.arg(tboot)
  modelo.orig <- x$x
  if (inherits(modelo.orig, "tvvar")) 
  {
    VAR <- modelo.orig
  } 
  else 
  {
    VAR <- eval.parent(modelo.orig)
  }
  p <- VAR$p
  neq <- VAR$neq
  X <- VAR$x
  obs <- VAR$obs
  type <- VAR$type
  B <- tvBcoef(VAR)
  BOOT <- vector("list", runs)
  ystar <- matrix(0, nrow = VAR$totobs, ncol = neq)
  y.names <- colnames(VAR$y)
  colnames(ystar) <- y.names
  Zdet <- NULL
  if (NCOL(X) > (neq * (p + 1)))
  {#in case there are exogen variables
    Zdet <- as.matrix(X[, (neq * (p + 1) + 1):NCOL(X)])
  }
  beta.pilot <- NULL
  if(tboot == "lbwb") 
  {
    l <- floor(obs^(1/3))
    model.pilot <- update(VAR, bw =  2*(VAR$bw)^(5/9)) 
    resorig <- scale(residuals(model.pilot), scale = FALSE)
    fitted <- model.pilot$fitted
    beta.pilot <- model.pilot$coefficients
  } 
  else 
  {
    resorig <- scale(VAR$residuals, scale = FALSE)
    prob <- c(0.7236067977499789360962267892318777740001678466796875,
              0.2763932022500210639037732107681222259998321533203125)
    binar <- c((1-sqrt(5))*0.5, (1+sqrt(5))*0.5)
    fitted <- VAR$fitted
  }

  yorig <- VAR$y.orig
  for (i in 1:runs)
  {
    if (tboot == "wild")
    {
      resid <- resorig*sample(binar, obs, replace = TRUE, prob = prob)
    }
    else if (tboot == "wild2")
    {
      resid  <-resorig*stats::rnorm(obs*neq)
    }
    else if (tboot == "lbwb") 
    {
      # Multivariate LBWB: block the entire (neq) vector
      obs.blocks <- ceiling(obs/l)
      resid <- matrix(0, nrow = obs, ncol = neq)
      
      for (k in 1:obs.blocks)
      {
        l.index <- max(1, k * l - 2 * l + 1)
        u.index <- min(obs - l + 1, k * l + 1)
        pi.i <- sample(l.index:u.index, 1)
        block.selected <- resorig[pi.i:(pi.i + l - 1), , drop = FALSE]
        fill.start <- (k - 1) * l + 1
        len_to_fill <- min(l, obs - fill.start + 1)
        if(len_to_fill > 0) 
        {
          resid[fill.start:(fill.start + len_to_fill - 1), ] <- 
            stats::rnorm(1) * block.selected[1:len_to_fill, , drop = FALSE]
        }
      }
    }
    lasty <- c(t(yorig[p:1, ]))
    ystar[c(1:p), ] <- yorig[c(1:p), ]
    for (j in 1:obs)
    {
      lasty <- lasty[1:(neq * p)]
      Z <- c(lasty, Zdet[j, ])
      if (type == "const")
        Z <- c(Z, 1)
      ystar[j + p, ] <- B[j,,]%*%Z + resid[j, ]
      lasty <- c(ystar[j + p, ], lasty)
    }
    VAR$y.orig <- ystar
    VAR$y <- ystar[-c(1:p),]
    VAR$residuals <- resid
    temp <- stats::embed(ystar, dimension = p + 1)[, -c(1:neq)]
    if(type == "const")
      temp <- cbind(temp, 1L)
    VAR$x[, 1:NCOL(temp)] <- temp
    varboot <- update(VAR)
    BOOT[[i]] <- .tvIRF(x = varboot, impulse = impulse, response = response, y.names = y.names,
                         n.ahead = n.ahead, ortho = ortho, cumulative = cumulative,
                         ortho.cov = ortho.cov, bw.cov = bw.cov, unit.shock = unit.shock)$irf
  }
  return(list(BOOT = BOOT, beta.pilot = beta.pilot))
}
