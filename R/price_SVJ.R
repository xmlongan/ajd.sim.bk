#' Pricing the European Call Option Under the SVJ Model
#'
#' @description
#' Pricing the European call option under the SVJ model using the
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
#'  underling asset price process.
#' @param mu_bar parameter \eqn{\bar{\mu}}, see Broadie-Kaya (2006).
#' @param sigma_s parameter \eqn{\sigma_s}, see Broadie-Kaya (2006).
#' @param true_price theoretical true price of the option, calculable by
#'  the Bates (1996) method.
#' @param log.inv `FALSE`, logging the Fourier inversion if `log.inv = TRUE`.
#'
#' @return vector of pricing error and consumed time (simulation).
#' @export
#'
#' @examples
#' S = 100; K = 100; v0 = 0.008836; k = 3.99; theta = 0.014; sigma = 0.27
#' rho = -0.79; r = 0.0319; tau = 5; lambda = 0.11; mu_bar = -0.12
#' sigma_s = 0.15; true_price = 20.1642
#' # price_SVJ(10000, S, K, v0, tau, r, k, theta, sigma, rho, lambda, mu_bar,
#' # sigma_s, true_price)
price_SVJ  <- function(N, S, K, v0, tau, r, k, theta, sigma, rho,
                       lambda, mu_bar, sigma_s, true_price, log.inv=F) {
  start.time = Sys.time()
  Y = ryield_SVJ(N, v0, tau, r, k, theta, sigma, rho,
                 lambda, mu_bar, sigma_s, log.inv)
  cprice_MC = exp(-r*tau) * mean(pmax(S*exp(Y)-K, 0))
  end.time = Sys.time()
  time.taken = end.time - start.time
  error = cprice_MC - true_price
  #
  return(c(error, as.numeric(time.taken)))
}
