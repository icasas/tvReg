#' Create List of Control Parameters for tvSURE
#'
#' Create a list of control pararameters for function \code{\link{tvSURE}}.
#' All control parameters that are not passed to this function are set to
#' default values.
#'
#' If the estimation is iterative FGLS with \code{maxiter}>1, the convergence criterion is
#' \deqn{\sqrt{ \frac{ \sum_{i, j}
#' (B_{i,j,g} - B{i, j, g-1})^2 }{ \sum_{i, j} B_{i, j, g-1}^2 }} < \code{tol}}
#' (\eqn{B_{i, j,g}} is the ith, jth coefficient of the gth iteration step).
#' @rdname tvReg-internals
#' @param maxiter maximum number of iterations for the iterative FGLS
#' estimations.
#' @param tol tolerance level indicating when to stop the iteration for the iterative FGLS estimations
#' @return A list of the above components.
#'
#' @export
tvsure.control <- function(maxiter = 1, tol = 1e-05)
{
  result <- list()
  if (maxiter <= 0 || round(maxiter) != maxiter) {
    stop("control parameter 'maxiter' must be a positive integer\n")
  }
  result$maxiter <- maxiter
  if (tol <= 0 || !is.numeric(tol) || length(tol) != 1) {
    stop("control parameter 'tol' must be a positive scalar\n")
  }
  result$tol <- tol
  return(result)
}
