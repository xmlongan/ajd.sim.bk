#' Characteristic Function of IV
#'
#' @description
#' The conditional Characteristic Function (CF) of the integrated variance,
#' given the two side levels of the variance.
#' - `CF()`: Compute the CF values all at once.
#' - `CF_next()`: Compute the CF value one by one.
#'
#' @details
#' The integrated variance over period \eqn{(u,t]} is defined as
#' \deqn{
#'   IV_{u,t} \triangleq  \int_{u}^{t}v(s)ds.
#' }
#' The conditional CF is defined as
#' \eqn{
#'   \Phi(a) \triangleq  E\left[ e^{iaIV_{u,t}}|v_u,v_t \right],
#' }
#' which can be decomposed into three parts:
#' \eqn{
#'   \Phi(a) = f_1(a) \times f_2(a) \times f_3(a),
#' }
#' where
#' \deqn{
#'   f_1(a) \triangleq
#'   \frac{\gamma(a)}{k} e^{-0.5(\gamma(a)-k)(t-u)}
#'   \frac{1-e^{-k(t-u)}}
#'      {1-e^{-\gamma(a)(t-u)}},
#' }
#' \deqn{
#'   f_2(a) \triangleq
#'    \text{exp}\left\{ \frac{v_u+v_t}{\sigma^2}\left[
#'     k\frac{1+e^{-k(t-u)}}{1-e^{-k(t-u)}} -
#'     \gamma(a)\frac{1+e^{-\gamma(a)(t-u)}}{1-e^{-\gamma(a)(t-u)}}
#'   \right] \right\},
#' }
#' \deqn{
#'   f_3(a) \triangleq
#'    I_{0.5d-1}\left[
#'      \sqrt{v_uv_t}\frac{4\gamma(a)}{\sigma^2} \frac{e^{-0.5\gamma(a)(t-u)}}
#'          {1-e^{-\gamma(a)(t-u)}} \right]
#'                      \Big/I_{0.5d-1}\left[
#'      \sqrt{v_uv_t}\frac{4k}{\sigma^2} \frac{e^{-0.5k(t-u)}}
#'          {1-e^{-k(t-u)}} \right],
#' }
#' \eqn{\gamma(a) = \sqrt{k^2 - 2\sigma^2ai}}, \eqn{d = 4\theta k/\sigma^2},
#' and \eqn{I_v(\cdot)} is the modified Bessel function of the first kind.
#' - The evaluation of \eqn{I_v(z)} with a complex number input is done with
#' [Bessel::BesselI()] or [Bessel::besselIasym()], the numerator in \eqn{f_3(a)}.
#' - The evaluation of \eqn{I_v(x)} with a real number input is done with
#' [Bessel::BesselI()] or [Bessel::besselIasym()], the denominator in \eqn{f_3(a)}.
#'
#' The modified Bessel function of the first kind is given as the following
#' power series:
#' \deqn{
#'   I_v(z) = \left(\frac{1}{2}z\right)^v \sum_{j=0}^{\infty}
#'   \frac{\left(\frac{1}{4}z^2\right)^j}{j!\Gamma(v+j+1)},
#' }
#' where \eqn{\Gamma(x)} is the gamma function and \eqn{z} is a complex number.
#'
#' @name CF
#'
#' @param a vector of real numbers, does't work for a single number,
#' should also work for complex numbers.
#' @param v0 value of the left end, \eqn{v_u}.
#' @param v1 value of the right end, \eqn{v_t}.
#' @param tau time difference, \eqn{t-u}.
#' @param k parameter \eqn{k}.
#' @param theta parameter \eqn{\theta}.
#' @param sigma parameter \eqn{\sigma}.
#' @param a1 current input number, scalar, complex number should work still.
#' @param z0 previous input for the embedded function `cBesselI_next()`.
#' @param branch0 branch of z0.
#'
#' @return
#' - `CF()`: vector of complex numbers.
#' - `CF_next()`: list of (CF_next(a1), z1, branch1).
#' @export
#'
#' @references
#' Broadie, M., & Kaya, Ö. (2006). Exact simulation of stochastic volatility
#' and other affine jump diffusion processes. Operations research,
#' 54(2), 217-231.
#'
#' @examples
#' # a = seq(0, 1500, 0.001)
#' a = seq(0, 1500, 100)
#' k = 6.21; theta = 0.019; sigma = 0.61; v0 = 0.010201; v1 = v0; tau = 1
#' f = CF(a, v0, v1, tau, k, theta, sigma)
#' par(mfrow = c(2,1))
#' plot(Re(f), Im(f), type="l", col="blue")
#' plot(a, Re(f), type="l", main="Real part of the CF"); abline(h=0, col="red")
CF <- function(a, v0, v1, tau, k, theta, sigma) {
  if (length(a) < 2) {stop("Please use CF_next() for scalar input!")}
  ga = sqrt(k^2 - 2*sigma^2*a*1i) # gamma(a)
  dk = exp(-k*tau); da = exp(-ga*tau)
  #
  f = (ga/k) * exp(-0.5*(ga-k)*tau) * ((1-dk)/(1-da))
  #
  t1 = (v0+v1)/sigma^2
  t2 = k * ((1+dk)/(1-dk)) - ga * ((1+da)/(1-da))
  f = f * exp(t1 * t2)
  #
  nu = 2*theta*k/sigma^2 - 1
  zs = sqrt(v0*v1) * (4*ga/sigma^2) * exp(-0.5*ga*tau)/(1-da)
  # print(sprintf("z = (%.5f, %.5f i)", Re(zs[1]), Im(zs[1])))
  x  = sqrt(v0*v1) * (4*k  /sigma^2) * exp(-0.5*k*tau) / (1-dk)
  f = f * (cBesselI(zs, nu) / cBesselI(x,  nu))
  # zi = zs; zo = cBesselI(zi, nu); plot(Re(zo), Im(zo), type="l")
  return(f)
}
