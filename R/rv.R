#' Simulate Next Variance
#'
#' @description
#' Simulate the next variance given the current variance level.
#'
#' @param v0 current variance level.
#' @param tau time difference, \eqn{t-u}.
#' @param k parameter \eqn{k}.
#' @param theta parameter \eqn{\theta}.
#' @param sigma parameter \eqn{\sigma}.
#'
#' @return next variance, scalar.
#' @export
#'
#' @examples
#' k = 6.21; theta = 0.019; sigma = 0.61; v0 = 0.010201; tau = 1
#' v1 = rv(v0, tau, k, theta, sigma)
rv <- function(v0, tau, k, theta, sigma) {
  d = 4*theta*k/sigma^2             # degrees of freedom
  C = sigma^2*(1-exp(-k*tau))/(4*k) # leading constant
  lambda = v0*exp(-k*tau)/C         # noncentrality parameter
  return(C * stats::rchisq(1, d, ncp = lambda))
}
