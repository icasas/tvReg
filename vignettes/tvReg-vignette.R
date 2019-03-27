## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(echo = TRUE, tidy = TRUE, fig.path = "fig/", fig.width = 6,
                      fig.height = 4, fig.align = "center")
library(tvReg)

## ----eval = TRUE---------------------------------------------------------
 ##Simulate a linear process with time-varying coefficient
 ##as functions of scaled time.
 library(tvReg)
 tau <- seq(1:1000)/1000
 beta <- data.frame(beta1 = sin(2 * pi * tau), beta2 = 2 * tau)
 X1 <- rnorm(1000)
 X2 <- rchisq(1000, df = 4)
 error <- rt(1000, df = 10)
 y <- apply(cbind(X1, X2) * beta, 1, sum) + error
 data <-data.frame(y = y, X1 = X1, X2 = X2)

 ##Estimate coefficients with lm and tvLM for comparison
 coef.lm <- stats::lm(y ~ 0 + X1 + X2, data = data)$coef
 model.tvlm <- tvLM(y ~ 0 + X1 + X2, data = data)

 ##Plot the estimates of beta1
 par (mar = c(4,4,1,1), oma = c(1,1,1,1))
 plot(tau,beta[, 1], type = "l", main = "", ylab = expression(beta[1]),
 xlab = expression(tau), ylim = range(beta[,1], model.tvlm$tvcoef[, 1]))
 abline(h = coef.lm[1], col = 2)
 lines(tau, model.tvlm$tvcoef[, 1], col=4)
 legend("topright", c(expression(beta[1]), "lm", "tvlm"), col = c(1, 2, 4), 
        bty = "n", lty = 1, cex = 0.8)

## ----eval = TRUE---------------------------------------------------------
 ##Data generation
 set.seed (42)
 z <- stats::arima.sim (n = 1000, list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)), sd = sqrt(0.1796))
 beta <- data.frame(beta1 = sin(2 * pi * z), beta2 = 2 * z)
 y <- apply(cbind(X1, X2) * beta, 1, sum) + error
 data<-data.frame(y = y, X1 = X1, X2 = X2, z = z)
 
 ##Coefficients estimation
 coef.lm <- stats::lm(y ~ 0 + X1 + X2, data = data)$coef
 model.tvlm2 <- tvLM(y ~ 0 + X1 + X2, z = z, data = data, bw = 0.4, est = "ll")
 
 ##Plotting the estimates of beta1
 sort.index <- sort.int(z, index.return = TRUE)$ix
 par (mar = c(4,4,1,1), oma = c(1,1,1,1))
 plot(z[sort.index], beta[sort.index, 1], type = "l", main = "",
 ylab = expression(beta[1]), xlab = expression(z[t]),
 ylim = range(beta[,1], model.tvlm2$tvcoef[, 1]))
 abline(h = coef.lm[1], col = 2)
 lines(z[sort.index], model.tvlm2$tvcoef[sort.index, 1], col = 4)
 legend("top", c(expression(beta[1]), "lm", "tvLM"), col = c(1, 2, 4), bty = "n", 
        lty = 1, cex = 0.5)

## ----eval = FALSE--------------------------------------------------------
#   ##Simulate an tvAR(2) process
#   tt<-(1:1000)/1000
#   beta <- cbind( 0.5 * cos (2 * pi * tt), (tt - 0.5)^2)
#   y<-numeric(1000)
#   y[1] <- 0.5
#   y[2] <- -0.2
#   ##y(t) = beta1(t) y(t-1) + beta2(t) y(t-2) + ut
#   for (t in 3:1000)
#   {
#     y[t] <- y[(t-1):(t-2)] %*% beta[t,] + rnorm(1)
#   }
#   Y <- tail (y, 500)
#  
#   ##Estimate coefficients of process Y with ar.ols and tvAR
#   ##and compare them in a plot
#   model.ar.2p <- ar.ols(Y, aic = FALSE, order = 2, intercept = FALSE, demean = FALSE)
#   model.tvAR.2p <- tvAR(Y, p = 2, type = "none", est = "ll")

## ----eval = FALSE--------------------------------------------------------
#  ##Simulate a AR(1) process with coefficients depending on z
#  z <- runif(2000, -1, 1)
#  beta <- (z - 0.5)^2
#  y <- numeric(2000)
#  y[1] <- 0.5
#  error <- rnorm(2000)
#   ##y(t) = beta1(z(t)) y(t-1) +  ut
#  for (t in 2:2000)
#  {
#   y[t] <- y[(t-1)] %*% beta[t] + error[t]
#  }
#  ##Remove initial conditions effects
#  Z <- z #tail (z, 1500)
#  Y <- y #tail (y, 1500)
#  #beta <- tail (beta, 1500)
#  ##Estimate coefficients of process Y with ar.ols and tvAR
#  ##and compare them in a plot
#  model.ar.1p <- ar.ols(Y, aic = FALSE, order = 1, intercept = FALSE, demean = FALSE)
#  model.tvAR.1p.z <- tvAR(Y, p = 1, z = Z, type = "none", est = "ll")

## ----eval = FALSE--------------------------------------------------------
#   ##Estimate only coefficient from odd lags and the intercept
#   tvAR.6p <- tvAR(y, p = 6, type = "const", fixed = c(NA, 0, NA, 0, NA, 0, NA), est = "ll")

## ------------------------------------------------------------------------
##Inflation rate, unemployment rate and treasury bill interest rate for the US,
##as used by Primiceri (2005).
data(usmacro, package = "bvarsv")
model.VAR <- vars::VAR(usmacro, p = 4, type = "const")
model.tvVAR <- tvVAR(usmacro, p = 4, type = "const")

## ------------------------------------------------------------------------
##Estimate a the impulse response functions with irf and tvIRF from 
##previous vector autoregressive models
model.irf <- vars::irf(model.VAR)
model.tvIRF <- tvIRF(model.tvVAR)

## ----eval = FALSE--------------------------------------------------------
#  ##Estimate a the tvIRF with time-varying covariance function
#  ##The estimation bandwidth might be automatically calculated or input
#  model.tvIRF2 <- tvIRF(model.tvVAR, cumulative = TRUE)

## ----eval = FALSE--------------------------------------------------------
#  data("Kmenta", package = "systemfit")
#  eqDemand <- consump ~ price + income
#  eqSupply <- consump ~ price + farmPrice + trend
#  system <- list(demand = eqDemand, supply = eqSupply)
#  
#  ##OLS estimation of a system
#  ols.fit <- systemfit::systemfit(system, method = "OLS", data = Kmenta)
#  
#  ##tvOLS estimation of a system with the local linear estimator
#  tvols.fit <- tvSURE(system, data = Kmenta,  est = "ll")

## ----eval = FALSE--------------------------------------------------------
#  ##FGLS estimation - SURE estimation
#  fgls1.fit <- systemfit::systemfit(system, data = Kmenta, method = "SUR")
#  
#  ##tvFGLS estimation - tvSURE estimation
#  tvfgls1.fit <- tvSURE(system, data = Kmenta, method = "tvFGLS")

## ----eval = FALSE--------------------------------------------------------
#  ##Iterative FGLS estimation - SUR estimation
#  fgls2.fit <- systemfit::systemfit(system, data = Kmenta, method = "SUR", maxit = 100)
#  
#  ##Iterative tvFGLS estimation - SURE estimation using the local linear
#  tvfgls2.fit <- tvSURE(system, data = Kmenta, method = "tvFGLS",
#  control = list(tol = 0.001, maxiter = 100))

## ----eval = FALSE--------------------------------------------------------
#  ##Estimation with 2 restrictions
#  Rrestr <- matrix(0, 2, 7)
#  Rrestr[1, 3] <-  1
#  Rrestr[1, 7] <- -1
#  Rrestr[2, 2] <- -1
#  Rrestr[2, 5] <-  1
#  qrestr <- c(0, 0.5)
#  tvfgls.rest <- tvSURE(system, data = Kmenta, method = "tvFGLS",
#  R = Rrestr, r = qrestr, bw = tvfgls1.fit$bw, bw.cov = tvfgls1.fit$bw.cov)

## ------------------------------------------------------------------------
library(MASS)
##Generate two independent (uncorrelated) series
y<-cbind(rnorm(200, sd= 4), rnorm(200, sd=1))
##Calculate the bandwidth
bw.cov<-bwCov(y)
##Estimate variance-variance matrix
Sigma.hat<- tvCov(y, bw=bw.cov)
##The first time point estimate
print(Sigma.hat[,,1])
##The mean over time of all estimates
print(apply(Sigma.hat, 1:2, mean))

##Generate two dependent variables with a covariance of -0.5
y <- mvrnorm(n = 200, mu = c(0,0), Sigma = cbind(c(1, -0.5), c(-0.5, 4)))
##Calculate the bandwidth
bw.cov <- bwCov(y)
##Estimation the variables variance-covariance matrix
Sigma.hat <- tvCov(y, bw = bw.cov)
##The first time point estimate
print(Sigma.hat[, , 1])

## ----eval = FALSE--------------------------------------------------------
#  ##Obtain the 90% confidence interval of the coefficients for an object of class tvlm
#  model.tvlm.90 <- confint (model.tvlm, level = 0.90, runs = 50)
#  
#  ##Obtain the 95% confidence interval of the same object. This will reused the resamples
#  ##of object model.tvlm done previously. So the second confidence interval calculation
#  ##is faster
#  model.tvlm.95 <- confint(model.tvlm.90)
#  
#  ##80% confidence interval using normal wild bootstrap for object of class tvar
#  ##with 200 bootstrap resamples
#  model.tvAR.ci <- confint(model.tvAR.1p.z, tboot = "wild2", level = 0.8, runs = 50)
#  
#  ##95% confidence interval using Mammen's wild bootstrap for object of class tvirf
#  model.tvIRF <- confint(model.tvIRF)
#  
#  ##50% confidence interval using Mammen's wild bootstrap for object of class tvsure
#  tvols.fit <- confint(tvols.fit, level = 0.5)

## ----eval = FALSE--------------------------------------------------------
#  ##Plot coefficient estimates and confidence intervals (if calculated) of objects of class tvlm
#  plot(model.tvlm.90)
#  
#  ##Plot coefficient estimates of objects of class tvar.
#  plot(model.tvAR.ci)
#  
#  ##Plot the fitted values and residuals of each equation in the model
#  plot(model.tvVAR)
#  
#  ##Plot the mean all tvIRF over time
#  plot(model.tvIRF)
#  
#  ##Plot the effect of a shock in the interest rates (tbi) on the inflation (inf)
#  ##at time 100
#  plot(model.tvIRF, obs.index = 100, impulse = "tbi", response = "inf")
#  
#  ##Plot the cumulative effect on a shock in short term interest rates (tbi)
#  ##on the inflation (inf)
#  plot <- plot(model.tvIRF2, impulse ="tbi", response = "inf")
#  
#  ##Plot coefficient estimates and confidence intervals (if calculated) of objects
#  ##of class tvsure
#  plot(tvols.fit)

