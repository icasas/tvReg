#' tvReg: Time-Varying Coefficient for Single and Multi-Equation Regressions
#' 
#' This package covers a large range of semiparametric regression methods with 
#' time-varying coefficients using nonparametric kernel smoothing for the estimation.
#' @section Main functions:
#' The five basic functions in this package are \code{\link{tvLM}}, \code{\link{tvAR}}, 
#' \code{\link{tvSURE}}, \code{\link{tvPLM}}, \code{\link{tvVAR}} and \code{\link{tvIRF}}. 
#' Moreover, this package provides the \code{\link{confint}}, \code{\link{fitted}}, 
#' \code{\link{forecast}}, \code{\link{plot}}, \code{\link{predict}}, \code{\link{print}}, 
#' \code{\link{resid}} and \code{\link{summary}} methods adapted to the class attributes 
#' of the \code{tvReg}. 
#' In addition, it includes bandwidth selection methods, time-varying variance-covariance 
#' estimators and four estimation procedures: the time-varying ordinary least squares, 
#' which are implemented in the \code{\link{tvOLS}} methods, the time-varying 
#' generalised least squares for a list of equations, which is implemented in the 
#' \code{\link{tvGLS}} methods, time-varying pooled and random effects estimators for 
#' panel data, which are implemented in the \code{\link{tvRE}} and the time-varying 
#' fixed effects estimator, which is implemente in the \code{\link{tvFE}}.
#' @section Further information:
#' Details on the theory and applications to finance and macroeconomics can be found 
#' in Casas and Fernandez-Casal (2019, \url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3363526}),
#' and in the package vignette \url{https://icasas.github.io/tvReg/articles/tvReg.html}.
#' @author Isabel Casas (\email{casasis@@gmail.com}),
#' Ruben Fernandez-Casal (\email{rubenfcasal@@gmail.com}).
#' @section Acknowledgments: 
#' Funded by the Horizon 2020. Framework Programme of the European Union.
#' @name tvReg-package
#' @aliases tvReg
#' @docType package
#' @import methods
#' @import Matrix
#' @import bvarsv
#' @importFrom stats confint
#' @importFrom MASS mvrnorm
#' @importFrom vars VAR irf
#' @importFrom systemfit systemfit
#' @importFrom plm pdata.frame pdim index
#' @references
#' Casas, I. and Fernandez-Casal, R., \emph{tvReg: Time-varying Coefficient Linear 
#' Regression for Single and Multi-Equations in R} (April 1, 2019). 
#' Available at SSRN: \url{https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3363526}.
NULL


#' Fama and French portfolio daily returns and factors for international markets.
#'
#' A dataset containing the returns of four portfolios ordered by size and book-to-market. The
#' four portfolios are SMALL/LoBM, SMALL/HiBM, BIG/LoBM and BIG/HiBM in four international
#' markets: North America (NA), Japan (JP), Asia Pacific (AP) and Europe (EU). It also contains the
#' Fama/French 5 factors for each of the markets.
#' @name FF5F
#' @docType data
#'
#' @references
#' Kennet R. French - Data Library (2017) 
#' http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html#International
#' 
#' Fama, E. and French, K. R (1993) Common risk factors in the returns on stocks and bonds,
#' \emph{Journal of Financial Economics}, 3-56.
#'
#' Fama, E. F. and French, K. R (2015) A five-factor asset pricing model,
#' \emph{Journal of Financial Economics}, 116, 1-22.
#'
#' @format A data frame with 314 rows and 41 variables.
#' \describe{
#'   \item{Date}{Date, months from July 1990 until August 2016}
#'   \item{NA.SMALL.LoBM}{Monthly returns of portfolio SMALL/LoBM in North American market}
#'   \item{NA.SMALL.HiBM}{Monthly returns of portfolio SMALL/HiBM in North American market}
#'   \item{NA.BIG.LoBM}{Monthly returns of portfolio BIG/LoBM in North American market}
#'   \item{NA.BIG.HiBM}{Monthly returns of portfolio BIG/HiBM in North American market}
#'   \item{NA.Mkt.RF}{North American market excess returns, i.e return of the 
#'   market - market risk free rate}
#'   \item{NA.SMB}{SMB (Small Minus Big) for the North American market}
#'   \item{NA.HML}{HML (High Minus Low) for the North American market}
#'   \item{NA.RMW}{RMW (Robust Minus Weak) for the North American market}
#'   \item{NA.CMA}{CMA (Conservative Minus Aggressive) for the North American market}
#'   \item{NA.RF}{North American risk free rate}
#'   \item{JP.SMALL.LoBM}{Monthly returns of portfolio SMALL/LoBM in Japanese market}
#'   \item{JP.SMALL.HiBM}{Monthly returns of portfolio SMALL/HiBM in Japanese market}
#'   \item{JP.BIG.LoBM}{Monthly returns of portfolio BIG/LoBM in Japanese market}
#'   \item{JP.BIG.HiBM}{Monthly returns of portfolio BIG/HiBM in Japanese market}
#'   \item{JP.Mkt.RF}{Japanese market excess returns, i.e return of the 
#'   market - market risk free rate}
#'   \item{JP.SMB}{SMB (Small Minus Big) for the Japanese market}
#'   \item{JP.HML}{HML (High Minus Low) for the Japanese market}
#'   \item{JP.RMW}{RMW (Robust Minus Weak) for the Japanese market}
#'   \item{JP.CMA}{CMA (Conservative Minus Aggressive) for the Japanese market}
#'   \item{JP.RF}{Japanese risk free rate}
#'   \item{AP.SMALL.LoBM}{Monthly returns of portfolio SMALL/LoBM in Asia Pacific market}
#'   \item{AP.SMALL.HiBM}{Monthly returns of portfolio SMALL/HiBM in Asia Pacific market}
#'   \item{AP.BIG.LoBM}{Monthly returns of portfolio BIG/LoBM in Asia Pacific market}
#'   \item{AP.BIG.HiBM}{Monthly returns of portfolio BIG/HiBM in Asia Pacific market}
#'   \item{AP.Mkt.RF}{Asia Pacific market excess returns, i.e return of the 
#'   market - maket risk free rate}
#'   \item{AP.SMB}{SMB (Small Minus Big) for the Asia Pacific market}
#'   \item{AP.HML}{HML (High Minus Low) for the Asia Pacific market}
#'   \item{AP.RMW}{RMW (Robust Minus Weak) for the Asia Pacific market}
#'   \item{AP.CMA}{CMA (Conservative Minus Aggressive) for the Asia Pacific market}
#'   \item{AP.RF}{Asia Pacific risk free rate}
#'   \item{EU.SMALL.LoBM}{Excess return of portfolio SMALL/LoBM in European market}
#'   \item{EU.SMALL.HiBM}{Excess return of portfolio SMALL/HiBM in European market}
#'   \item{EU.BIG.LoBM}{Excess return of portfolio BIG/LoBM in European market}
#'   \item{EU.BIG.HiBM}{Excess return of portfolio BIG/HiBM in European market}
#'   \item{EU.Mkt.RF}{European market excess returns, i.e returns of the 
#'   market - market risk free rate}
#'   \item{EU.SMB}{SMB (Small Minus Big) for the European market}
#'   \item{EU.HML}{HML (High Minus Low) for the European market}
#'   \item{EU.RMW}{RMW (Robust Minus Weak) for the European market}
#'   \item{EU.CMA}{CMA (Conservative Minus Aggressive) for the European market}
#'   \item{EU.RF}{European risk free rate}
#' }
#' @source \url{http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html}
#' @keywords datasets
NULL

#' Standarised rates from a currency portfolio.
#' 
#' Aslanidis and Casas (2013) consider a portfolio of daily US dollar exchange rates of the Australian
#' dollar (AUS), Swiss franc (CHF), euro (EUR), British pound (GBP), South African rand (RAND),
#' Brazilian real (REALB), and Japanese yen (YEN) over the period from January 1, 1999
#' until May 7, 2010 (T = 2856 observations). This dataset contains the standarised rates after
#' "devolatilisation", i.e. standarising the rates using a GARCH(1,1) estimate of the volatility.
#
#' @name CEES
#' @docType data
#'
#' @references
#' Aslanidis, N. and Casas, I. (2013) Nonparametric correlation models for portfolio allocation,
#' \emph{Journal of Banking \& Finance}, 37, 2268 - 2283.
#'
#'
#' @format A data frame with 2855 rows and 8 variables. Below the standarised rates of daily US dollar
#' exchange rates of
#' \describe{
#'  \item{Date}{Daily data from Jan 6, 1999 until May 7, 2010 - without
#'   weekends and days off}
#'   \item{AUS}{Australian dollar}
#'   \item{CHF}{Swiss franc}
#'   \item{EUR}{Euro}
#'   \item{GBP}{British pound}
#'   \item{RAND}{South African rand}
#'   \item{REALB}{Brazilian real}
#'   \item{YEN}{Japanese yen}
#' }
#' @keywords datasets
NULL


#' Daily realized variance
#' 
#' A dataset containing the daily realized variance, and some of its lags,
#' obtained from 1-minute close prices of the S\&P 500. Similar data has
#' been used in the HAR model in Corsi (2009), the HARQ and SHARQ models in
#' Bollerslev et al (2016) and the TVHARQ and TVSHARQ models in 
#' Casas et al (2018). The time period runs from Feb 1990 until 
#' Dec 2006.
#
#' @name RV
#' @docType data
#'
#' @references
#' 
#' Bollerslev, T., Patton, A. J. and Quaedvlieg, R. (2016) Exploiting the 
#' errors: A simple approach for improved volatility forecasting. 
#' \emph{Journal of Econometrics}, 192, 1-18.
#' 
#' Bollerslev, T., Tauchen, G. and Zhou, H. (2009) Expected stock returns 
#' and variance risk premia. \emph{The Review of Financial Studies}, 22, 44-63.
#' 
#' Casas, I., Mao, X. and Vega, H. (2018) Reexamining financial and economic 
#' predictability with new estimators of realized variance and variance 
#' risk premium. Url= http://pure.au.dk/portal/files/123066669/rp18_10.pdf
#' 
#' Corsi, F. (2009) A simple approximate long-memory model of realized 
#' volatility. \emph{Journal of Financial Econometrics}, 7, 174-196.
#'
#'
#' @format A data frame with 4264 rows and 6 variables.
#' \describe{
#'   \item{Date}{Daily data from Feb 1, 1990 until Dec 29, 2006 - without
#'   weekends and days off}
#'   \item{RV}{Daily realized variance at time t}
#'   \item{RV_lag}{Daily realized variance at time t-1}
#'   \item{RV_week}{Weekly average realized variance at time t-5}
#'   \item{RV_month}{Monthly average realized variance at time t-22}
#'   \item{RQ_lag_sqrt}{Daily squared root of the realized quarticity at time t-1}
#' }
#' @keywords datasets
NULL

#' Variables related to the problem of healthcare spending. 
#' 
#' @name OECD
#' @docType data
#'
#' @references
#' Casas, I., Gao, J., Peng B., and Xie, S. (2019). Ferreira, E., and Orbe, S. (2017) 
#' Modelling Time-Varying Income Elasticities of Health Care Expenditure for the OECD. 
#' Available at SSRN: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3262326
#' 
#' @format A data frame with 680 rows and 7 columns. 
#' \describe{
#'   \item{country}{}
#'   \item{year}{}
#'   \item{lhe}{Log of country's healthcare spending}
#'   \item{lgdp}{log of country's gdp}
#'   \item{pop65}{Country's ratio of population greater than 65 years old}
#'   \item{pop14}{Country's ratio of population younger than 15 years old}
#'   \item{public}{Country's ratio of healthcare funding coming from the government}
#' }
#' @keywords datasets
NULL