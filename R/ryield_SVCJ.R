#' Simulate Yield using Broadie-Kaya Method for the SVCJ model
#'
#' @description
#' Simulate the next period yield given current variance using the Broadie-Kaya
#' method for the SVCJ model.
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
#' underling asset price process and variance process.
#' @param mu_bar parameter \eqn{\bar{\mu}}.
#' @param sigma_s parameter \eqn{\sigma_s}.
#' @param mu_v parameter \eqn{\mu_v}, mean of the jumps
#' (exponential distribution) in the variance.
#' @param rho_J parameter \eqn{\rho_J}, correlation parameter between the
#' jumps in the price and variance.
#' @param log.inv `FALSE`, logging the Fourier inversion if `log.inv = TRUE`.
#'
#' @return next period yield, a scalar.
#' @export
#'
#' @examples
#' v0 = 0.007569; k = 3.46; theta = 0.008; sigma = 0.14; rho = -0.82
#' r = 0.0319; tau = 1; lambda = 0.47; mu_bar = -0.1; sigma_s = 0.0001
#' mu_v = 0.05; rho_J = -0.38
#' # Y = ryield_Hest(1000, v0, tau, r, k, theta, sigma, rho,
#' # lambda, mu_bar, sigma_s, mu_v, rho_J)
#' # hist(Y)
ryield_SVCJ <- function(n, v0, tau, r, k, theta, sigma, rho,
                        lambda, mu_bar, sigma_s, mu_v, rho_J,
                        log.inv = FALSE) {
  # change the r (mu)
  r = r - lambda * mu_bar
  mu_s = log((1+mu_bar)*(1-rho_J*mu_v)) - sigma_s^2/2
  #
  Y = rep(0, n)
  #
  v0_original = v0
  #
  if (log.inv) {
    fname = format(Sys.time(), "./script/logs/riv-%Hhour-%Mmin.csv")
    if (!file.exists(fname)) {
      title = c("v0", "v1", "tau", "k", "theta", "sigma", "h",
                "num_Bessel_eval", "secs_consumed\n")
      cat(title, file=fname, sep=',', append=T)
    }
  }
  #
  for (i in 1:n) {
    v0 = v0_original
    #
    numJ = stats::rpois(1, lambda * tau)
    # no jump
    if (numJ == 0) {
      Y[i] = ryield_Hest(1, v0, tau, r, k, theta, sigma, rho)
      next
    }
    # at least one jump
    Jtime    = stats::runif(numJ, 0, tau)
    delta_ts = diff(c(0, sort(Jtime))) # Jtime unchanged after sorting
    for (j in 1:numJ) {
      # diffusion
      delta_t = delta_ts[j]
      v1 = rv(v0, delta_t, k, theta, sigma)
      iv = riv(1, v0, v1, delta_t, k, theta, sigma, log.inv)
      I  = (v1 - v0 - k*theta*delta_t + k*iv)/sigma
      mu = r*delta_t - iv/2 + rho*I
      diffusion = mu + sqrt(1-rho^2) * stats::rnorm(1, mean=0, sd=sqrt(iv))
      # jump
      J_v = stats::rexp(1, rate = 1/mu_v)
      v1 = v1 + J_v
      # J_s = exp(stats::rnorm(1, mean=mu_s+rho_J*J_v, sd=sigma_s)) - 1
      J_s = stats::rnorm(1, mean=mu_s+rho_J*J_v, sd=sigma_s)
      #
      Y[i] = Y[i] + diffusion + J_s
      #
      v0 = v1
    }
    # residual time diffusion
    delta_t = tau - Jtime[numJ]
    if (delta_t > 0) {
      Y[i] = Y[i] + ryield_Hest(1, v0, delta_t, r, k, theta, sigma, rho)
    }
  }
  return(Y)
}
