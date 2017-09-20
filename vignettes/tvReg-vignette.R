## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, tidy = TRUE)
library(tvReg)

## ------------------------------------------------------------------------
library( "tvReg" )
## Data simulation
tau <- seq(1:1000)/1000
beta <- data.frame(beta1 = sin(2 * pi * tau), beta2 = 2 * tau)
X1 <- rnorm(1000)
X2 <- rchisq(1000, df = 4)
error <- rt(1000, df = 10)
y <- apply(cbind(X1, X2) * beta, 1, sum) + error
data<-data.frame(y = y, X1 = X1, X2 = X2)

## Estimate coefficients with lm and tvLM for comparison

coef.lm <- stats::lm(y~0+X1+X2, data = data)$coef
coef.tvlm <- tvLM(y~0+X1+X2, data = data)$tvcoef

## Comparing results in a plot
plot(tau,beta[, 1], type = "l", main = "", ylab = expression(beta[1]),
xlab = expression(tau), ylim = range(beta[,1], coef.tvlm[, 1]))
abline(h = coef.lm[1], col = 2)
lines(tau, coef.tvlm[, 1], col = 4)
legend("topright", c(expression(beta[1]), "lm", "tvlm"),
col = c(1, 2, 4), bty = "n", lty = 1)

## ------------------------------------------------------------------------
library(MASS)
##Generate two independent (uncorrelated) series
y<-cbind(rnorm(200, sd= 4), rnorm(200, sd=1))
## Calculate the bandwidth
bw.cov<-bwCov(y)
## Estimate variance-variance matrix
Sigma.hat<- tvCov(y, bw=bw.cov)
## The first time point estimate
print(Sigma.hat[,,1])
## The mean over time of all estimates
print(apply(Sigma.hat, 1:2, mean))

## ------------------------------------------------------------------------
## Generate two dependent variables with a covariance of -0.5
y <- mvrnorm(n = 200, mu = c(0,0), Sigma = cbind(c(1, -0.5), c(-0.5, 4)))
## Calculate the bandwidth
bw.cov <- bwCov(y)
## Estimation the variables variance-covariance matrix
Sigma.hat <- tvCov(y, bw = bw.cov)
## The first time point estimate
print(Sigma.hat[, , 1])

