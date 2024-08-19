#' Truncate the Characteristic Function
#'
#' @description
#' Truncate the Characteristic Function automatically through setting up the
#' truncation error bound.
#'
#' @param h grid size.
#' @param v0 the left side variance level, \eqn{v_u}.
#' @param v1 the right side variance level, \eqn{v_t}.
#' @param tau time difference, \eqn{t-u}.
#' @param k parameter \eqn{k}.
#' @param theta parameter \eqn{\theta}.
#' @param sigma parameter \eqn{\sigma}.
#' @param eps truncate error bound.
#'
#' @return
#' vector of CF points.
#' @export
#'
#' @examples
#' k = 6.21; theta = 0.019; sigma = 0.61; v0 = 0.010201; v1 = v0; tau = 1
#' mu = miv(v0, v1, tau, k, theta, sigma); sd = sqrt(mu[2] - mu[1]^2)
#' mu_eps = mu[1] + 7*sd; h = pi/mu_eps
#' zs = CF_trunc(h, v0, v1, tau, k, theta, sigma)
#' plot(Re(zs), Im(zs), type="b", col="blue"); abline(h=0,v=0, col="red")
#' hj = seq(1, length(zs)) * h
#' plot(hj, Re(zs), type="b", col="blue"); abline(h=0, col="red")
CF_trunc <- function(h, v0, v1, tau, k, theta, sigma, eps=1e-5) {
  # dk = exp(-k*tau)
  # z0 = sqrt(v0*v1)*(4*k/sigma^2)*exp(-0.5*k*tau)/(1-dk); branch0 = 0
  #
  z0 = 1+1i; branch0 = 0 # any value in the first quadrant is OK for z0
  f_z1_branch1 = CF_next(h, z0, branch0, v0, v1, tau, k, theta, sigma)
  f = f_z1_branch1$f; z1 = f_z1_branch1$z1; branch1 = f_z1_branch1$branch1
  #
  j = 1
  fs = c(f)
  #
  while (Mod(f)/j >= pi*eps/2) {
    j = j+1; z0 = z1; branch0 = branch1
    f_z1_branch1 = CF_next(j*h, z0, branch0, v0, v1, tau, k, theta, sigma)
    f = f_z1_branch1$f; z1 = f_z1_branch1$z1; branch1 = f_z1_branch1$branch1
    fs = append(fs, f)
  }
  return(fs)
}
