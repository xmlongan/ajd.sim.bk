#' Conditional Distribution of IV
#'
#' @description
#' Density, distribution function and random variable generation for the
#' conditional distribution of the Integrated Variance (IV),
#' all through Fourier inversion of the Characteristic Function,
#' see [CF()] for the details.
#'
#' Currently, only [riv()] has been tested heavily and can be used
#' safely for most scenarios.
#'
#' While, for other functions,
#' use them carefully, because currently the truncation of Fourier inversion of
#' the CF is not done in a robust way.
#'
#' @param q vector of quantile.
#' @param x vector of quantile.
#' @param n number of observations.
#' @param v0 the left end, \eqn{v_u}.
#' @param v1 the right end, \eqn{v_t}.
#' @param tau time difference, \eqn{t-u}.
#' @param k parameter \eqn{k}.
#' @param theta parameter \eqn{\theta}.
#' @param sigma parameter \eqn{\sigma}.
#' @param log.inv `FALSE`, logging the Fourier inversion if `log.inv = TRUE`.
#'
#' @details
#' The distribution function is inversed from the conditional CF, see
#' [invCF()] and [CF()].
#' Noting that `piv()` is a wrapper for `invCF()`.
#' Currently, the discretized inversions are done through:
#' - `piv()`, `div()`:  `hj = seq(0.001, 1500, 0.001)`, i.e., `h=0.001`
#' and truncated at `j=length(hj)`.
#' - `riv()`: \eqn{h = \pi/\mu_{\epsilon}} with
#' \eqn{\mu_{\epsilon} = \mu_{iv} + 7\sigma_{iv}}, where \eqn{\mu_{iv}} and
#' \eqn{\sigma_{iv}} denote the conditional mean and standard deviation of the
#' integrated variance, respectively. Meanwhile, the truncation is done at
#' \eqn{j=N} when
#' \deqn{
#'   \frac{|\Phi(hN)|}{N} < \frac{\pi\epsilon}{2},
#' }
#' where \eqn{\epsilon} is the desired truncation error.
#'
#' @return
#' `div` gives the density, `piv` gives the distribution function, and
#' `riv` generates random deviates.
#' @export
#'
#' @seealso See [CF()] for the details about the Characteristic
#' Function.
#'
#' @references
#' Broadie, M., & Kaya, Ö. (2006). Exact simulation of stochastic volatility
#' and other affine jump diffusion processes. Operations research,
#' 54(2), 217-231.
#'
#' @examples
#' # x = seq(0,0.1,0.001)
#' # k = 6.21; theta = 0.019; sigma = 0.61; v0 = 0.010201; v1 = v0; tau = 1
#'
#' # Pr = piv(x, v0, v1, tau, k, theta, sigma)
#' # par(mfrow=c(1,1))
#' # txt = "Conditional CDF of IV"
#' # plot(x,Pr,type="l",col="blue",ylab="CDF",main=txt)
#' # abline(h=0,v=0); abline(h=1, col="red")
#'
#' # d = div(x, v0, v1, tau, k, theta, sigma)
#' # txt = "Conditional PDF of IV"
#' # plot(x, d, type="l", col="blue", ylab="PDF", main=txt); abline(h=0,v=0)
#'
#' # r = riv(1, v0, v1, tau, k, theta, sigma)
#'
piv <- function(q, v0, v1, tau, k, theta, sigma) {
  return(invCF(q, v0, v1, tau, k, theta, sigma))
}
