#' Pricing the European Call Option Under the Heston Model
#'
#' @description
#' Pricing the European call option under the Heston model using
#' Broadie-Kaya (2006) method for the simulation.
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
#' @param true_price theoretical price of the option, calculable from
#' the Heston (1993) method.
#' @param log.inv `FALSE`, logging the Fourier inversion if `log.inv = TRUE`.
#'
#' @return vector of pricing error and consumed time (simulation).
#' @export
#'
#' @examples
#' S = 100; K = 100; v0 = 0.010201; k = 6.21; theta = 0.019; sigma = 0.61
#' rho = -0.7; r = 0.0319; tau = 1; true_price = 6.8061
#' # price_Hest(10000, S, K, v0, tau, r, k, theta, sigma, rho, true_price)
price_Hest <- function(N, S, K, v0, tau, r, k, theta, sigma, rho,
                       true_price, log.inv = F) {
  start.time = Sys.time()
  Y = ryield_Hest(N, v0, tau, r, k, theta, sigma, rho, log.inv)
  cprice_MC = exp(-r*tau) * mean(pmax(S*exp(Y)-K, 0))
  end.time = Sys.time()
  time.taken = end.time - start.time
  error = cprice_MC - true_price
  #
  return(c(error, as.numeric(time.taken)))
}
