#' Continued BesselI function
#'
#' @description
#' Accomplish the continuity of the BesselI (the modified Bessel function
#' of the first kind) function. see [CF()].
#' The following continuation formula from Abramowitz and Stegun (1972):
#' \deqn{
#'   I_v(ze^{m\pi i}) = e^{mv\pi i} I_v(z), \{m~~\text{integer}\},
#' }
#' is used to calculate the value of \eqn{I_v(z)} for \eqn{arg(z)} different
#' than its principle value.
#'
#' @param zs vector of complex or real numbers.
#' @param v numeric (scalar), i.e., the order (maybe fractional and
#' negative) of the Modified Bessel function of the first kind, \eqn{v} in
#' \eqn{I_v(z)}.
#' @param z1 current input complex number.
#' @param z0 previous input complex number.
#' @param branch0 branch for z0.
#' @param expon_scaled whether the BesselI() should be scaled.
#'
#' @return
#' - `cBesselI()`: vector of complex or real numbers.
#' - `cBesselI_next()`: list of current input's function value and branch.
#' @export
#'
#' @examples
#' # a = seq(0, 1500, 0.001)
#' a = seq(0, 1500, 100)
#' k = 6.21; theta = 0.019; sigma = 0.61; v0 = 0.010201; v1 = v0; tau = 1
#' ga = sqrt(k^2 - 2*sigma^2*a*1i)
#' dk = exp(-k*tau); da = exp(-ga*tau)
#' zs = sqrt(v0*v1) * (4*ga/sigma^2) * exp(-0.5*ga*tau)/(1-da)
#' nu = 2*theta*k/sigma^2 - 1
#' f = cBesselI(zs, nu)
#' par(mfrow=c(2,1))
#' plot(Re(zs), Im(zs), col="blue", type="l"); abline(h=0, v=0, col="red")
#' plot(Re(f), Im(f), col="blue", type="l", main="Continued BesselI function")
#'
#' z1 = -1 - 0.1i; z0 = -1 + 0.1i; branch0 = 0
#' f_branch = cBesselI_next(z1, z0, branch0, nu)
cBesselI <- function(zs, v) {
  branches = get_branch(zs); N = length(branches); f = rep(0,N)
  for (n in 1:N) {
    z = zs[n]; branch = branches[n]
    if (is.complex(z)) {
      cond1 = abs(Re(z)) < 1000; cond2 = abs(Re(z)) > 400
    } else {
      cond1 = abs(z) < 1000; cond2 = abs(z) > 400
    }
    mybesselI = ifelse(cond1, Bessel::BesselI, Bessel::besselIasym)
    expon_scaled = ifelse(cond2, TRUE, FALSE)
    f[n] = exp(v*branch*1i) * mybesselI(z, v, expon.scaled = expon_scaled)
  }
  return(f)
}
