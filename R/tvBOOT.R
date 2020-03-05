# @rdname tvReg-internals
#' @keywords internal
#'
.tvboot<-function (x , ...) UseMethod(".tvboot", x)

# @rdname tvReg-internals
#' @method .tvboot default
#' @keywords internal
#'
.tvboot.default<-function (x, runs = 100, tboot = "wild")
{
  if (!(any(class(x) %in% c("tvlm", "tvar"))))
  {
    stop("Bootstrap confidence intervals not implemented for this class.\n")
  }
  obs <- length(x$y)
  nvar <- NCOL (x$x)
  resorig <- scale(x$residuals, scale = FALSE)
  fitted <- x$fitted
  residup <- resorig * -0.6180339887498949025257 # (1-sqrt(5))*0.5
  residdown <- resorig * 1.618033988749894902526 # (1+sqrt(5))*0.5
  yboot <- matrix(NA, nrow = obs, ncol = runs)
  for (i in 1:runs) 
  {
    if (tboot == "wild"){
      index <- ifelse(apply(resorig, 2, stats::rbinom, size = 1, prob = 0.7236067977499789360962) == 1, TRUE, FALSE)
      resid <- residup
      resid [!index] <- residdown[!index]
      ystar <- fitted + resid
    }else if (tboot == "wild2"){
      resid <- resorig*stats::rnorm(obs)
      ystar <- fitted + resid
    }
    yboot[, i] <- ystar
  }
  return(.tvLM.ci(x, yboot))
}

# @rdname tvReg-internals
#' @method .tvboot tvplm
#' @keywords internal
#'
.tvboot.tvplm <- function (x , runs = 100, tboot = "wild")
{
  B <- x$coefficients
  BOOT <- vector("list", runs)
  resorig <- matrix(x$residuals, nrow = x$obs, ncol = x$neq)
  resorig = scale(resorig, scale = FALSE)
  residup <- resorig * -0.6180339887498949025257 # (1-sqrt(5))*0.5
  residdown <- resorig * 1.618033988749894902526 # (1+sqrt(5))*0.5
  fitted <- x$fitted
  X.tilde <- x$x
  bw <- x$bw
  obs <- x$obs
  neq <- x$neq
  yboot <- matrix(NA, nrow = obs*neq, ncol = runs)
  for (i in 1:runs){
    if (tboot=="wild"){
      index <- ifelse(apply(resorig, 2, stats::rbinom, size=1, prob=0.7236067977499789360962)==1, TRUE, FALSE)
      resid <- residup
      resid [!index] <- residdown[!index]
      ystar <- fitted + as.numeric(resid)
    }else if (tboot=="wild2"){
      resid <- resorig*stats::rnorm(obs*neq)
      ystar <- fitted + resid
    }
    yboot[, i] <- ystar
  }
  return(.tvPLM.ci(x, yboot))
}

# @rdname tvReg-internals
#' @method .tvboot tvsure
#' @keywords internal
#'
.tvboot.tvsure<-function (x, runs = 100, tboot = "wild")
{
  obs <- x$obs
  nvar<- x$nvar
  neq <- x$neq
  resorig <- scale(x$residuals, scale = FALSE)
  fitted <- x$fitted
  residup <- resorig * -0.6180339887498949025257 # (1-sqrt(5))*0.5
  residdown <- resorig * 1.618033988749894902526 # (1+sqrt(5))*0.5
  yboot <- vector("list", runs)
  for (i in 1:runs) {
    if (tboot == "wild"){
      index <- ifelse(apply(resorig, 2, stats::rbinom, size = 1, prob = 0.7236067977499789360962) == 1, TRUE, FALSE)
      resid <- residup
      resid [!index] <- residdown[!index]
      ystar <- fitted + resid
    }else if (tboot == "wild2"){
      resid <- resorig*stats::rnorm(obs*neq)
      ystar <- fitted + resid
    }
    yboot[[i]] <- ystar
  }
  return(.tvSURE.ci(x, yboot))
}

# @rdname tvReg-internals
#' @method .tvboot tvirf
#' @keywords internal
.tvboot.tvirf<-function (x, runs = 0, tboot = "wild")
{
  ortho <- x$ortho
  cumulative <- x$cumulative
  impulse <- x$impulse
  response <- x$response
  bw.cov <- x$bw.cov
  ortho.cov <- x$ortho.cov
  n.ahead <- x$n.ahead
  x <- x$x
  if (inherits(x, "tvvar"))
  {
    VAR <- eval.parent(x)
  }
  else
  {
    stop("Bootstrap not implemented for this class.\n")
  }
  p <- VAR$p
  neq <- VAR$neq
  X <- VAR$x
  obs <- VAR$obs
  total <- VAR$totobs
  type <- VAR$type
  params <- NCOL(x$x)
  B <- tvBcoef(VAR)
  BOOT <- vector("list", runs)
  ystar <- matrix(0, nrow = total, ncol = neq)
  y.names <- colnames(VAR$y)
  colnames(ystar) <- y.names
  Zdet <- NULL
  if (NCOL(X) > (neq * (p + 1)))
  {#in case there are exogen variables
    Zdet <- as.matrix(X[, (neq * (p + 1) + 1):NCOL(X)])
  }
  resorig <- scale(VAR$residuals, scale = FALSE)
  Sigmaorig <- crossprod(resorig)/(obs - params)
  eSigma <- eigen(Sigmaorig, TRUE)
  Sigma <- diag(eSigma$values^0.5)%*%eSigma$vector
  resid <-resorig
  fitted <- VAR$fitted
  yorig <- VAR$y.orig
  residup <- resorig * -0.6180339887498949025257 # (1-sqrt(5))*0.5
  residdown <- resorig * 1.618033988749894902526 # (1+sqrt(5))*0.5
  for (i in 1:runs)
  {
    if (tboot == "wild")
    {
      index <- ifelse(apply(resorig, 2, stats::rbinom, size=1, prob=0.7236067977499789360962) == 1, TRUE, FALSE)
      resid <- residup
      resid [!index] <- residdown[!index]
      lasty <- c(t(yorig[p:1, ]))
      ystar[c(1:p), ] <- yorig[c(1:p), ]
    }
    else if (tboot == "wild2")
    {
      resid  <-resorig*stats::rnorm(obs*neq)
      lasty <- c(t(yorig[p:1, ]))
      ystar[c(1:p), ] <- yorig[c(1:p), ]
    }
    for (j in 1:obs)
    {
      lasty <- lasty[1:(neq * p)]
      Z <- c(lasty, Zdet[j, ])
      if (VAR$type == "const")
        Z <- c(Z, 1)
      ystar[j + p, ] <- B[j,,]%*%Z + resid[j, ]
      lasty <- c(ystar[j + p, ], lasty)
    }
    VAR$y.orig <- ystar
    VAR$y <- ystar[-c(1:p),]
    VAR$residuals <- resid
    temp <- stats::embed(ystar, dimension = p + 1)[, -c(1:neq)]
    if(VAR$type == "const")
      temp <- cbind(temp, rep(1, obs))
    VAR$x[, 1:NCOL(temp)] <- temp
    varboot <- update(VAR)
    BOOT[[i]] <- .tvIRF(x = varboot, impulse = impulse, response = response, y.names = y.names,
                         n.ahead = n.ahead, ortho = ortho, cumulative = cumulative,
                         ortho.cov = ortho.cov, bw.cov = bw.cov)$irf
  }
  return(BOOT)
}
