#' Cumulative Distribution Function of IV
#'
#' @description
#' Get the conditional Cumulative Distribution Function (CDF) of the
#' integrated variance, given the two side variance levels.
#' Use it carefully, because currently the truncation of Fourier inversion of
#' the CF is not done in a robust way.
#'
#' @details
#' The CDF is computed using Fourier inversion methods:
#' \deqn{
#'   F(x)\equiv \text{P}\{IV_{u,t} \le x\} = \frac{2}{\pi}\int_{0}^{\infty}
#'   \frac{\sin(ux)}{u}\text{Re}[\Phi(u)]du.
#' }
#' The trapezoidal rule is used to compute the probability distribution
#' numerically:
#' \deqn{
#'   \text{P}\{IV_{u,t} \le x\} \approx \frac{hx}{\pi} + \frac{2}{\pi}
#'   \sum_{j=1}^{\infty}\frac{\sin(hjx)}{j}\text{Re}[\Phi(hj)],
#' }
#' where \eqn{h} is the grid size, which is set as 0.001 internally, and
#' `hj = seq(0.001, 1500, 0.001)`.
#'
#' @param x a scalar or vector of real numbers, support of the CDF.
#' @param v0 the left end, \eqn{v_u}.
#' @param v1 the right end, \eqn{v_t}.
#' @param tau time difference, \eqn{t-u}.
#' @param k parameter \eqn{k}.
#' @param theta parameter \eqn{\theta}.
#' @param sigma parameter \eqn{\sigma}.
#'
#' @return probability/probabilities corresponding to \eqn{x}, depending on
#' whether the supplied \eqn{x} is a scalar or a vector.
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
#' k = 6.21; theta = 0.019; sigma = 0.61; v0 = 0.010201; v1 = v0; tau = 1
#' # Pr = invCF(x, v0, v1, tau, k, theta, sigma)
#' # par(mfrow=c(1,1))
#' # txt = "Cumulative Probability Function"
#' # plot(x, Pr, type="l", col="blue", ylab="CDF", main=txt)
#' # abline(h=0, v=0); abline(h=1, col="red")
invCF <- function(x, v0, v1, tau, k, theta, sigma) {
  # need to first plot the curve of CF, to get a feeling of the variation
  hj = seq(0.001, 1500, 0.001); j = 1:length(hj); h = hj[2] - hj[1]
  # sensible to the initial value
  zs = CF(hj, v0, v1, tau, k, theta, sigma)
  #
  N = length(x)
  if (N == 1) {
    prob = h*x/pi + (2/pi) * sum((sin(hj*x)/j) * Re(zs))
  } else {
    prob = rep(0, N)
    for (i in 1:N) {
      prob[i] = h*x[i]/pi + (2/pi) * sum((sin(hj*x[i])/j) * Re(zs))
    }
  }
  return(prob)
}
