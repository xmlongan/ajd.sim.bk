#' @rdname piv
#' @export
div <- function(x, v0, v1, tau, k, theta, sigma) {
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
  N = length(x); d = rep(0, N)
  for (i in 1:N) {
    #     Set up grid size to bound discretization error
    h = 2*pi/(x[i] + mu_eps)
    #      Truncate the Characteristic Function inversion
    zs = CF_trunc(h, v0, v1, tau, k, theta, sigma, eps=1e-5)
    j = seq(1, length(zs)); hj = h*j; Rezs = Re(zs); rm(zs)
    #
    d[i] = h/pi + (2/pi) * sum(cos(hj*x[i]) * h * Re(zs))
  }
  return(d)
}
