#' Simulate Yield using Broadie-Kaya Method for the SVJ model
#'
#' @description
#' Simulate the next period yield given current variance using
#' Broadie-Kaya (2006) method for the SVJ model.
#'
#' @param n number of yield samples to simulate.
#' @param v0 current variance level.
#' @param tau time difference, \eqn{t-u}.
#' @param r riskless rate, risk-neutral drift.
#' @param k parameter \eqn{k}.
#' @param theta parameter \eqn{\theta}.
#' @param sigma parameter \eqn{\sigma}.
#' @param rho parameter \eqn{\rho}.
#' @param lambda parameter \eqn{\lambda}, arrival rate of jumps in the
#'   underling asset price process.
#' @param mu_bar parameter \eqn{\bar{\mu}}.
#' @param sigma_s parameter \eqn{\sigma_s}.
#' @param log.inv `FALSE`, logging the Fourier inversion if `log.inv = TRUE`.
#'
#' @return next period yield, a scalar.
#' @export
#'
#' @examples
#' v0 = 0.008836; k = 3.99; theta = 0.014; sigma = 0.27; rho = -0.79
#' r = 0.0319; tau = 5; lambda = 0.11; mu_bar = -0.12; sigma_s = 0.15
#' # Y = ryield_SVJ(1000, v0, tau, r, k, theta, sigma, rho,
#' # lambda, mu_bar, sigma_s)
#' # hist(Y)
ryield_SVJ <- function(n, v0, tau, r, k, theta, sigma, rho,
                       lambda, mu_bar, sigma_s, log.inv = FALSE) {
  # diffusion part
  Y = ryield_Hest(n, v0, tau, r-lambda*mu_bar, k, theta, sigma, rho, log.inv)
  # jump part
  Jsum = rep(0, n)
  mu_s = log(1+mu_bar) - sigma_s^2/2
  Ns = stats::rpois(n, lambda*tau)
  for (i in 1:n) {
    if (Ns[i] > 0) {
      # Jsum[i] = sum(exp(stats::rnorm(Ns[i], mu_s, sigma_s)) - 1)
      Jsum[i] = sum(stats::rnorm(Ns[i], mu_s, sigma_s))
    }
  }
  # continuous + jump
  return(Y+Jsum)
}
