[![CRAN Status](https://www.r-pkg.org/badges/version/tvReg)](https://cran.r-project.org/package=tvReg)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/last-month/tvReg?color=orange)](https://cran.r-project.org/package=tvReg)


# tvReg

This R package covers a large range of semiparametric regression methods with time-varying coefficients using nonparametric kernel smoothing for the estimation.


## Installation

You can install the released version of tvReg from [CRAN](https://CRAN.R-project.org) with:

```
install.packages("tvReg")
```

or the development version from GitHub with:
```
devtools::install_github("icasas/tvReg")
```

## Main functions

The five basic functions in this package are `tvLM()`, `tvAR()`, `tvSURE()`, `tvPLM()`, `tvVAR()` and `tvIRF()`. Moreover, this package provides the `confint()`, `fitted()`, `forecast()`, `plot()`, `predict()`, `print()`, `resid()` and `summary()` methods adapted to the class attributes of the `tvReg`. In addition, it includes bandwidth selection methods, time-varying variance-covariance estimators and four estimation procedures: the time-varying ordinary least squares, which are implemented in the `tvOLS()` methods, the time-varying generalised least squares for a list of equations, which is implemented in the `tvGLS()` methods, time-varying pooled and random effects estimators for panel data, which are implemented in the `tvRE()` and the time-varying fixed effects estimator, which is implemente in the `tvFE()`.


## References

Casas, Isabel and Fernandez-Casal, Ruben, *tvReg: Time-varying Coefficient Linear Regression for Single and Multi-Equations in R* (2022). R Journal, 14/1, pp. 79 - 100. [link](https://journal.r-project.org/articles/RJ-2022-002/).



