#' Fama and French portfolio excess returns and factors for international markets.
#'
#' A dataset containing the excess returns of four portfolios ordered by size and book-to-market. The
#' four portfolios are SMALL/LoBM, SMALL/HiBM, BIG/LoBM and BIG/HiBM in four international
#' market: North America (NA), Japan (JP), Asia Pacific (AS) and Europe (EU). It also contains the
#' Fama/French 5 factors for each of the markets.
#' @name FF5F
#' @docType data
#'
#' @references
#' Fama, E. and French, K. R (1993) Common risk factors in the returns on stocks and bonds,
#' \emph{Journal of Financial Economics}, 3-56.
#'
#' Fama, E. F. and French, K. R (2015) A five-factor asset pricing model,
#' \emph{Journal of Financial Economics}, 116, 1-22.
#'
#' @format A data frame with 314 rows and 41 variables.
#' \describe{
#'   \item{Date}{Date, months from July 1990 until Agust 2016}
#'   \item{NA.SMALL.LoBM}{Excess return of portfolio SMALL/LoBM in North American market}
#'   \item{NA.SMALL.HiBM}{Excess return of portfolio SMALL/HiBM in North American market}
#'   \item{NA.BIG.LoBM}{Excess return of portfolio BIG/LoBM in North American market}
#'   \item{NA.BIG.HiBM}{Excess return of portfolio BIG/HiBM in North American market}
#'   \item{NA.Mkt.RF}{North American market excess return}
#'   \item{NA.SMB}{SMB (Small Minus Big) for the North American market}
#'   \item{NA.HML}{HML (High Minus Low) for the North American market}
#'   \item{NA.RMW}{RMW (Robust Minus Weak) for the North American market}
#'   \item{NA.CMA}{CMA (Conservative Minus Aggressive) for the North American market}
#'   \item{NA.RF}{North American risk free rate}
#'   \item{JP.SMALL.LoBM}{Excess return of portfolio SMALL/LoBM in Japanese market}
#'   \item{JP.SMALL.HiBM}{Excess return of portfolio SMALL/HiBM in Japanese market}
#'   \item{JP.BIG.LoBM}{Excess return of portfolio BIG/LoBM in Japanese market}
#'   \item{JP.BIG.HiBM}{Excess return of portfolio BIG/HiBM in Japanese market}
#'   \item{JP.Mkt.RF}{Japanese market excess return}
#'   \item{JP.SMB}{SMB (Small Minus Big) for the Japanese market}
#'   \item{JP.HML}{HML (High Minus Low) for the Japanese market}
#'   \item{JP.RMW}{RMW (Robust Minus Weak) for the Japanese market}
#'   \item{JP.CMA}{CMA (Conservative Minus Aggressive) for the Japanese market}
#'   \item{JP.RF}{Japanese risk free rate}
#'   \item{AP.SMALL.LoBM}{Excess return of portfolio SMALL/LoBM in Asia Pacific market}
#'   \item{AP.SMALL.HiBM}{Excess return of portfolio SMALL/HiBM in Asia Pacific market}
#'   \item{AP.BIG.LoBM}{Excess return of portfolio BIG/LoBM in Asia Pacific market}
#'   \item{AP.BIG.HiBM}{Excess return of portfolio BIG/HiBM in Asia Pacific market}
#'   \item{AP.Mkt.RF}{Asia Pacific market excess return}
#'   \item{AP.SMB}{SMB (Small Minus Big) for the Asia Pacific market}
#'   \item{AP.HML}{HML (High Minus Low) for the Asia Pacific market}
#'   \item{AP.RMW}{RMW (Robust Minus Weak) for the Asia Pacific market}
#'   \item{AP.CMA}{CMA (Conservative Minus Aggressive) for the Asia Pacific market}
#'   \item{AP.RF}{Asia Pacific risk free rate}
#'   \item{EU.SMALL.LoBM}{Excess return of portfolio SMALL/LoBM in European market}
#'   \item{EU.SMALL.HiBM}{Excess return of portfolio SMALL/HiBM in European market}
#'   \item{EU.BIG.LoBM}{Excess return of portfolio BIG/LoBM in European market}
#'   \item{EU.BIG.HiBM}{Excess return of portfolio BIG/HiBM in European market}
#'   \item{EU.Mkt.RF}{European market excess return}
#'   \item{EU.SMB}{SMB (Small Minus Big) for the European market}
#'   \item{EU.HML}{HML (High Minus Low) for the European market}
#'   \item{EU.RMW}{RMW (Robust Minus Weak) for the European market}
#'   \item{EU.CMA}{CMA (Conservative Minus Aggressive) for the European market}
#'   \item{EU.RF}{European risk free rate}
#' }
#' @source \url{http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html}
#' @keywords datasets
NULL

#' Aslanidis and Casas (2013) consider a portfolio of daily US dollar exchange rates of the Australian
#' dollar (AUS), Swiss franc (CHF), euro (EUR), British pound (GBP), South African rand (RAND),
#' Brazilian real (REALB), and Japanese yen (YEN) over the period from January 1, 1999
#' until May 7, 2010 (T=2856 observations). This dataset contains the standarised rates after
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
#' @format A data frame with 2856 rows and 7 variables. Below the standarised rates of daily US dollar
#' exchange rates of
#' \describe{
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
