#' Pricing the European Call Option Under the SVCJ Model
#'
#' @description
#' Pricing the European call option under the SVCJ model using the
#' Broadie-Kaya method for the simulation.
#'
#' @param N number of underling asset price samples to simulate.
#' @param S current underling asset price.
#' @param K striking price.
#' @param v0 current variance level.
#' @param tau maturity time.
#' @param r riskless rate.
#' @param k parameter \eqn{k}.
#' @param theta parameter \eqn{\theta}.
#' @param sigma parameter \eqn{\sigma}.
#' @param rho parameter \eqn{\rho}.
#' @param lambda parameter \eqn{\lambda}, arrival rate of jumps in the
#' underling asset price process and variance process.
#' @param mu_bar parameter \eqn{\bar{\mu}}, see Broadie-Kaya (2006).
#' @param sigma_s parameter \eqn{\sigma_s}, see Broadie-Kaya (2006).
#' @param mu_v parameter \eqn{\mu_v}, mean of the jumps
#' (exponential distribution) in the variance.
#' @param rho_J parameter \eqn{\rho_J}, correlation parameter between the
#' jumps in the price and variance.
#' @param true_price theoretical true price of the option, calculable by
#' the Duffie et al. (2000) method.
#' @param log.inv `FALSE`, logging the Fourier inversion if `log.inv = TRUE`.
#'
#' @return vector of pricing error and consumed time (simulation).
#' @export
#'
#' @examples
#' S = 100; K = 100; v0 = 0.007569; k = 3.46; theta = 0.008; sigma = 0.14
#' rho = -0.82; r = 0.0319; tau = 1; lambda = 0.47; mu_bar = -0.1
#' sigma_s = 0.0001; mu_v = 0.05; rho_J = -0.38; true_price = 6.8619
#' # price_SVCJ(10000, S, K, v0, tau, r, k, theta, sigma, rho, lambda, mu_bar,
#' # sigma_s, mu_v, rho_J, true_price)
price_SVCJ <- function(N, S, K, v0, tau, r, k, theta, sigma, rho,
                       lambda, mu_bar, sigma_s, mu_v, rho_J, true_price,
                       log.inv=F) {
  start.time = Sys.time()
  Y = ryield_SVCJ(N, v0, tau, r, k, theta, sigma, rho,
                  lambda, mu_bar, sigma_s, mu_v, rho_J, log.inv)
  cprice_MC = exp(-r*tau) * mean(pmax(S*exp(Y)-K, 0))
  end.time = Sys.time()
  time.taken = end.time - start.time
  error = cprice_MC - true_price
  #
  return(c(error, as.numeric(time.taken)))
}
