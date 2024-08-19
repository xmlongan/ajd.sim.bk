#' Derivative of the Density of IV
#'
#' @description
#' Derivative of the conditional Density of the Integrated Variance (IV).
#' Use it carefully, because currently the truncation of Fourier inversion of
#' the CF is not done in a robust way.
#'
#' @param x vector of quantile.
#' @param v0 the left end, \eqn{v_u}.
#' @param v1 the right end, \eqn{v_t}.
#' @param tau time difference, \eqn{t-u}.
#' @param k parameter \eqn{k}.
#' @param theta parameter \eqn{\theta}.
#' @param sigma parameter \eqn{\sigma}.
#'
#' @return
#' derivative/derivatives corresponding to \eqn{x}, depending on
#' whether the supplied \eqn{x} is a scalar or a vector.
#' @export
#'
#' @examples
#' # x = seq(0,0.1,0.001)
#' # k = 6.21; theta = 0.019; sigma = 0.61; v0 = 0.010201; v1 = v0; tau = 1
#'
#' # dd = ddiv(x, v0, v1, tau, k, theta, sigma)
#' # txt = "Derivative of Conditional PDF of IV"
#' # plot(x, dd, type="l", col="blue", ylab="Derivative", main=txt)
#' # abline(h=0,v=0)
ddiv <- function(x, v0, v1, tau, k, theta, sigma) {
  hj = seq(0.001, 1500, 0.001); j = 1:length(hj); h = hj[2] - hj[1]
  zs = CF(hj, v0, v1, tau, k, theta, sigma)
  #
  N = length(x)
  if (N == 1) {
    dd = - (2/pi) * sum(sin(hj*x) * hj * h * Re(zs))
  } else {
    dd = rep(0, N)
    for (i in 1:N) {
      dd[i] = - (2/pi) * sum(sin(hj*x[i]) * hj * h * Re(zs))
    }
  }
  return(dd)
}
