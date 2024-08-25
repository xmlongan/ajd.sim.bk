#' Simulate Yield using Broadie-Kaya Method for the Heston Model
#'
#' @description
#' Simulate the next period yield given current variance using
#' Broadie-Kaya (2006) method.
#'
#' @param n number of yield samples to simulate.
#' @param v0 current variance level.
#' @param tau time difference, \eqn{t-u}.
#' @param r riskless rate, risk-neutral drift.
#' @param k parameter \eqn{k}.
#' @param theta parameter \eqn{\theta}.
#' @param sigma parameter \eqn{\sigma}.
#' @param rho parameter \eqn{\rho}.
#' @param log.inv `FALSE`, logging the Fourier inversion if `log.inv = TRUE`.
#'
#' @return next period yield, a scalar.
#' @export
#'
#' @examples
#' v0 = 0.010201; k = 6.21; theta = 0.019; sigma = 0.61; rho = -0.7
#' r = 0.0319; tau = 1
#' # Y = ryield_Hest(1000, v0, tau, r, k, theta, sigma, rho)
#' # hist(Y)
ryield_Hest <- function(n, v0, tau, r, k, theta, sigma, rho, log.inv = F) {
  #
  Y = rep(0, n)
  for (i in 1:n) {
    v1 = rv(v0, tau, k, theta, sigma)
    iv = riv(v0, v1, tau, k, theta, sigma, log.inv)
    I = (v1 - v0 - k*theta*tau + k*iv)/sigma
    mu = r*tau - iv/2 + rho*I
    Y[i] = mu + sqrt(1-rho^2) * stats::rnorm(1, 0, sd=sqrt(iv))
  }
  return(Y)
}
