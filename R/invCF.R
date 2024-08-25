#' Cumulative Distribution Function of IV
#'
#' @description
#' Get the conditional Cumulative Distribution Function (CDF) of the
#' integrated variance, given the two side variance levels.
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
#' where \eqn{h} is the grid size, set through
#' \eqn{h = 2\pi/(x + \mu_{\epsilon}}) with
#' \eqn{\mu_{\epsilon} = \mu_{iv} + 5\sigma_{iv}}, where \eqn{\mu_{iv}} and
#' \eqn{\sigma_{iv}} denote the conditional mean and standard deviation of the
#' integrated variance, respectively. Meanwhile, the truncation is done at
#' \eqn{j=N} when
#' \deqn{
#'   \frac{|\Phi(hN)|}{N} < \frac{\pi\epsilon}{2},
#' }
#' where \eqn{\epsilon} is the desired truncation error.
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
  #                       Cautions
  # when tau is close to 0, it will produce large x or z
  # for the modified Bessel function of the first kind I_v(x).
  # Consequently, approximation errors would probably occur.
  #                   Special Preprocessing
  # when tau is very small, assume IV is fixed
  if (tau < 0.0001) {
    stop("Numerical error would happen, so we don't implement for small tau.")
  }
  #            Compute the first two moments
  mu = miv(v0, v1, tau, k, theta, sigma)
  if (any(!is.finite(mu))) {
    stop(sprintf("mu[1]=%.20f, mu[2]=%.20f",mu[1],mu[2]))
  }
  var = mu[2] - (mu[1])^2
  if (var <= 0) {
    temp = "miv(%.10f, %.10f, %.10f, %.4f, %.4f, %.4f)"
    msg1 = sprintf(temp, v0, v1, tau, k, theta, sigma)
    temp = "var(iv) = %.20f, mu[1]=%.20f, mu[2]=%.20f"
    msg2 = sprintf(temp, var, mu[1], mu[2])
    cat(paste("var(iv) <= 0", msg1, msg2, "", sep="\n"))
    stop("Numerical error happens for the supplied parameters")
  }
  #              discretization error setting
  sd = sqrt(var); mu_eps = mu[1] + 5*sd
  #
  N = length(x); prob = rep(0, N)
  for (i in 1:N) {
    #     Set up grid size to bound discretization error
    h = 2*pi/(x[i] + mu_eps)
    #      Truncate the Characteristic Function inversion
    zs = CF_trunc(h, v0, v1, tau, k, theta, sigma, eps=1e-5)
    j = seq(1, length(zs)); hj = h*j; Rezs = Re(zs); rm(zs)
    #
    prob[i] = h*x[i]/pi + (2/pi) * sum((sin(hj*x[i])/j) * Re(zs))
  }
  return(prob)
}
