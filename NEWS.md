# tvReg 0.3.0

* Add class `tvlm` to ojects of class `tvar`.

* `confint.tvlm()`, `confint.tvirf()` and `confint.tvsure()` added for compatibility with R standards

* `CI.tvlm()`, `CI.tvar()`, `CI.tvsure()` and `CI.tvirf()` deprecated

* `coef.tvlm()`, `coef.tvsure()`, `coef.tvvar()` and `coef.tvirf()` added for compatility with R standards

* `resid.tvlm()`, `resid.tvsure()`, `resid.tvvar()` and `resid.tvirf()` added for compatility with R standards

* `print.tvlm()`, `print.tvsure()`, `print.tvvar()` and `print.tvirf()` added for compatility with R standards

* `summary.tvlm()`, `summary.tvsure()`, `summary.tvvar()` and `summary.tvirf()` added for compatility with R standards

* `plot.tvsure()`, `plot.tvirf()` and `plot.tvlm()` to have an homegeneous tvReg plot format.

* Solve bug in `tvSURE()`. The use of parameter `est = "ll"` did not have effect. Previously, it always run the function as if `est = "lc"`.


