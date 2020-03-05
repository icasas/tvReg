#tvReg 0.5.0

* Added methods tvPLM, tvFE and tvRE methods to fit time-varying coefficients panel data models. Also added their corresponding confint, predict, forecast, plot, summary and print methods.

* Added dataset OECD

* Change all return argument "tvcoef" to "coefficients" to fit the standards of other packages in R

* Allow the user to choose individual plots for each variable in models of class tvlm, tvar or tvplm

* Fix of the bug in R-devel caused by matrix objects now also inheriting from class "array" which caused some problems in some "if" statements

#tvReg 0.4.2

* Added  modified or leave-(2l+1)-out cross-validation for bandwidth selection

* Faster algorithm for the calculation of confidence intervals

* Fix bug in `.tvCov.cv

# tvReg 0.4.1

* Fix bug in `tvGLS.R` for option `est = "ll"`.

* Added a website for the package (with pkgdown).

* Added `README.Rmd` and `index.Rmd`.

* Updated documentation `tvOLS matrix`, `forecast.tvar`, `forecast.tvlm` and
  `forecast.tvsure`.
  

# tvReg 0.4.0

* Removed class `tvlm` from objects of class `tvar`.

* Removed deprecated `CI.tvlm()`, `CI.tvar()`, `CI.tvsure()` and `CI.tvirf()` methods.

* Added `predict.tvlm()`, `predict.tvar()`, `predict.tvsure()` and `predict.tvvar()` methods.

* Added `forecast.tvlm()`, `forecast.tvar()`, `forecast.tvsure()` and `forecast.tvvar()` methods.

* Added `resid.tvlm()`, `resid.tvsure()`, `resid.tvvar()` and `resid.tvirf()` for compatility with R standards.



