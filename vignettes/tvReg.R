## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 6,
                      fig.height = 4, fig.align = "center")
library(tvReg)

## ----tvLM---------------------------------------------------------------------
 ##Simulate a linear process with time-varying coefficient
 ##as functions of rescaled time period.
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
 model.tvLM <- tvLM(y ~ 0 + X1 + X2, data = data)

 ##Plot the estimates of beta1
 par (mar = c(4,4,1,1), oma = c(1,1,1,1))
 plot(tau,beta[, 1], type = "l", main = "", ylab = expression(beta[1]),
 xlab = expression(tau), ylim = range(beta[,1], model.tvLM$coefficients[, 1]))
 abline(h = coef.lm[1], col = 2)
 lines(tau, model.tvLM$coefficients[, 1], col=4)
 legend("topright", c(expression(beta[1]), "lm", "tvlm"), col = c(1, 2, 4), 
        bty = "n", lty = 1, cex = 0.8)

## ----CI, eval = TRUE----------------------------------------------------------

##Obtain the 90% confidence interval of the coefficients for an object of the class attribute tvlm
model.tvLM.90 <- confint (model.tvLM, level = 0.90, runs = 50)

##Obtain the 95% confidence interval of the same object. This will reused the resamples of object model.tvLM.90. So the second confidence interval calculation is faster
model.tvLM.95 <- confint(model.tvLM.90)

## ----PLOT,  eval = TRUE-------------------------------------------------------
##Plot coefficient estimates and confidence intervals (if calculated) of objects of the class attribute tvlm
plot(model.tvLM.90)

## ----tvLM2--------------------------------------------------------------------
 ##Data generation
 set.seed (42)
 z <- stats::arima.sim (n = 1000, list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)), sd = sqrt(0.1796))
 beta <- data.frame(beta1 = sin(2 * pi * z), beta2 = 2 * z)
 y <- apply(cbind(X1, X2) * beta, 1, sum) + error
 data<-data.frame(y = y, X1 = X1, X2 = X2, z = z)
 
 ##Coefficients estimation
 coef.lm <- stats::lm(y ~ 0 + X1 + X2, data = data)$coef
 model.tvLM2 <- tvLM(y ~ 0 + X1 + X2, z = z, data = data, cv.block = 50, est = "ll")
 
 ##Plotting the estimates of beta1
 sort.index <- sort.int(z, index.return = TRUE)$ix
 par (mar = c(4,4,1,1), oma = c(1,1,1,1))
 plot(z[sort.index], beta[sort.index, 1], type = "l", main = "",
 ylab = expression(beta[1]), xlab = expression(z[t]),
 ylim = range(beta[,1], coefficients(model.tvLM2)[, 1]))
 abline(h = coef.lm[1], col = 2)
 lines(z[sort.index], model.tvLM2$coefficients[sort.index, 1], col = 4)
 legend("top", c(expression(beta[1]), "lm", "tvLM"), col = c(1, 2, 4), bty = "n", 
        lty = 1, cex = 0.5)

## ----tvAR---------------------------------------------------------------------
 ##Simulate an tvAR(2) process
 tt<-(1:1000)/1000
 beta <- cbind( 0.5 * cos (2 * pi * tt), (tt - 0.5)^2)
 y<-numeric(1000)
 y[1] <- 0.5
 y[2] <- -0.2
 ##y(t) = beta1(t) y(t-1) + beta2(t) y(t-2) + ut
 for (t in 3:1000)
 {
   y[t] <- y[(t-1):(t-2)] %*% beta[t,] + rnorm(1)
 }
 Y <- tail (y, 500)
 
 ##Coefficient estimates of process Y with ar.ols and tvAR
 model.ar.2p <- ar.ols(Y, aic = FALSE, order = 2, intercept = FALSE, demean = FALSE)
 model.tvAR.2p <- tvAR(Y, p = 2, type = "none", est = "ll")

## ----tvAR.z-------------------------------------------------------------------
##Simulate a AR(1) process with coefficients depending on z
z <- runif(2000, -1, 1)
beta <- (z - 0.5)^2
y <- numeric(2000)
y[1] <- 0.5
error <- rnorm(2000)
 ##y(t) = beta1(z(t)) y(t-1) +  ut
for (t in 2:2000)
{
 y[t] <- y[(t-1)] %*% beta[t] + error[t]
}

##Remove initial conditions effects
Z <- tail (z, 1500)
Y <- tail (y, 1500)

##Coefficient estimates of process Y with ar.ols and tvAR
model.ar.1p <- ar.ols(Y, aic = FALSE, order = 1, intercept = FALSE, demean = FALSE)
model.tvAR.1p.z <- tvAR(Y, p = 1, z = Z, type = "none", est = "ll")

##80% confidence interval using normal wild bootstrap for object of the class attribute tvar with 200 bootstrap resamples 
model.tvAR.80 <- confint(model.tvAR.1p.z, tboot = "wild2", level = 0.8, runs = 50)

##Plot coefficient estimates of objects of the class attribute tvar.
plot(model.tvAR.80)


## ----PRINT, eval = TRUE-------------------------------------------------------
##Summary of model.tvAR.80
summary(model.tvAR.80)

##Print of model.tvAR.80
print(model.tvAR.80)

## ----tvAR2, eval = FALSE------------------------------------------------------
#   ##Estimate only coefficient from odd lags and the intercept
#   tvAR.6p <- tvAR(y, p = 6, type = "const", fixed = c(NA, NA, 0, NA, 0, NA, 0), est = "ll")

## ----tvOLS--------------------------------------------------------------------
data("Kmenta", package = "systemfit")
eqDemand <- consump ~ price + income
eqSupply <- consump ~ price + farmPrice + trend
system <- list(demand = eqDemand, supply = eqSupply)

##OLS estimation of a system
model.ols <- systemfit::systemfit(system, method = "OLS", data = Kmenta)

##tvOLS estimation of a system with the local linear estimator
model.tvOLS <- tvSURE(system, data = Kmenta,  est = "ll")

##50% confidence interval using Mammen's wild bootstrap for object of the class attribute tvsure
model.tvOLS <- confint(model.tvOLS, level = 0.5, runs = 30)

##Plot coefficient estimates and confidence intervals (if calculated) of objects of the class attribute tvsure
plot(model.tvOLS)

##Print of model model.tvOLS 
print(model.tvOLS)

## ----tvFGLS-------------------------------------------------------------------
##FGLS estimation - SURE estimation
fgls1.fit <- systemfit::systemfit(system, data = Kmenta, method = "SUR")

##tvFGLS estimation - tvSURE estimation
tvfgls1.fit <- tvSURE(system, data = Kmenta, method = "tvFGLS")

## ----tvFGLS2------------------------------------------------------------------
##Iterative FGLS estimation - SUR estimation
fgls2.fit <- systemfit::systemfit(system, data = Kmenta, method = "SUR", maxit = 100)

##Iterative tvFGLS estimation - SURE estimation using the local linear
tvfgls2.fit <- tvSURE(system, data = Kmenta, method = "tvFGLS",
control = list(tol = 0.001, maxiter = 100))

## ----tvFGLS.rest,  eval = FALSE-----------------------------------------------
#  ##Estimation with 2 restrictions
#  Rrestr <- matrix(0, 2, 7)
#  Rrestr[1, 7] <- 1
#  Rrestr[2, 2] <- 1
#  Rrestr[2, 5] <- -1
#  qrestr <- c(0, 0.5)
#  tvfgls.rest <- tvSURE(system, data = Kmenta, method = "tvFGLS",
#  R = Rrestr, r = qrestr, bw = tvfgls1.fit$bw, bw.cov = tvfgls1.fit$bw.cov)

## ----tvPLM,  eval = TRUE------------------------------------------------------
data("OECD")
print(names(OECD))
elast.fe <- plm::plm(lhe ~ lgdp + pop65 + pop14 + public, data = OECD, 
                     index = c("country", "year"), model = "within")
elast.tvfe <- tvPLM (lhe ~ lgdp + pop65 + pop14 + public, data = OECD, 
                     index = c("country", "year"), method = "within", bw = 0.6)
elast.fe <- confint(elast.fe)
elast.tvfe <- confint(elast.tvfe)
plot(elast.tvfe, vars = 1, ylim = c(0.45, 1.3))
graphics::par(mfrow = c(1, 1), 
                mar = c(4, 4, 2, 1), oma = c(0, 0, 0, 0))
x.axis <- 1:elast.tvfe$obs
par(new = TRUE)
plot(x.axis, rep(mean(elast.fe[1,]), elast.tvfe$obs), ylim = c(0.45, 1.3),
     ylab ="", xlab ="", type = "l")
graphics::polygon(c(rev(x.axis), x.axis), 
                  c(rep(rev(elast.fe[1, 2]), elast.tvfe$obs),
                    rep(elast.fe[1, 1], elast.tvfe$obs)),
                  col = "grey10", border = NA, fillOddEven = TRUE)
lines(x.axis, rep(mean(elast.fe[1,]), elast.tvfe$obs), col ="white")

## ----tvVAR--------------------------------------------------------------------
##Inflation rate, unemployment rate and treasury bill interest rate for the US,
##as used by Primiceri (2005).
data(usmacro, package = "bvarsv")
model.VAR <- vars::VAR(usmacro, p = 4, type = "const")
model.tvVAR <- tvVAR(usmacro, p = 4, type = "const")

##Plot the fitted values and residuals of each equation in the model 
plot(model.tvVAR)
 

## ----tvIRF--------------------------------------------------------------------
##Estimate a the impulse response functions with irf and tvIRF from 
##previous vector autoregressive models
model.irf <- vars::irf(model.VAR)
model.tvIRF <- tvIRF(model.tvVAR)

##95% confidence interval using Mammen's wild bootstrap for object of the class attribute tvirf
model.tvIRF <- confint(model.tvIRF, runs = 30)

## ----tvIRF2-------------------------------------------------------------------
##Plot the mean all tvIRF over time
plot(model.tvIRF)

##Plot the effect of a shock in the interest rates (tbi) on the inflation (inf)
##at time 100
plot(model.tvIRF, obs.index = 100, impulse = "tbi", response = "inf")

## ----tvIRF3-------------------------------------------------------------------
##Estimate a the tvIRF with time-varying covariance function
model.tvIRF2 <- tvIRF(model.tvVAR, cumulative = TRUE)

##Plot the cumulative effect on a shock in short term interest rates (tbi) 
##on the inflation (inf)
plot <- plot(model.tvIRF2, impulse ="tbi", response = "inf")



## ----Sigma--------------------------------------------------------------------
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

## ----Forecast, eval = TRUE----------------------------------------------------
data(RV)
RV2 <- head(RV, 2001)
##Estimate/train tvHAR model
tvHAR <- with(RV2, tvAR (RV, p = 1, bw = 0.8, exogen = cbind(RV_week, RV_month)))
##Define the forecast horizon (n.ahead) and the future values of the exogenous variables
newexogen <- cbind(RV$RV_week[2002:2004], RV$RV_month[2002:2004])
##3-step-ahead forecast
forecast(tvHAR, n.ahead = 3, newexogen = newexogen)

## ----Predict, eval = TRUE-----------------------------------------------------
tvHARQ <- with(RV2, tvLM (RV ~ RV_lag + RV_week + RV_month, z = RQ_lag_sqrt, 
                bw = 0.003))
newdata <- cbind(RV$RV_lag[2002:2004], RV$RV_week[2002:2004], RV$RV_month[2002:2004])
newz <- RV$RQ_lag_sqrt[2002:2004]
predict(tvHARQ, newdata, newz)

