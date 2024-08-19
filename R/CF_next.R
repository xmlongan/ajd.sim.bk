#' @rdname CF
#' @export
CF_next <- function(a1, z0, branch0, v0, v1, tau, k, theta, sigma) {
  a = a1; if (length(a) > 1) {stop("Please use CF() for vector input!")}
  ga = sqrt(k^2 - 2*sigma^2*a*1i) # gamma(a)
  # print(sprintf("ga = (%.5f, %.5fi)", Re(ga), Im(ga)))
  dk = exp(-k*tau); da = exp(-ga*tau)
  #
  f = (ga/k) * exp(-0.5*(ga-k)*tau) * ((1-dk)/(1-da))
  #
  t1 = (v0+v1)/sigma^2
  t2 = k * ((1+dk)/(1-dk)) - ga * ((1+da)/(1-da))
  f = f * exp(t1 * t2)
  #
  nu = 2*theta*k/sigma^2 - 1
  #
  z1 = sqrt(v0*v1) * (4*ga/sigma^2) * exp(-0.5*ga*tau)/(1-da)
  x  = sqrt(v0*v1) * (4*k  /sigma^2) * exp(-0.5*k*tau) / (1-dk)
  #
  num_scaled = abs(Re(z1)) > 400
  val_branch = cBesselI_next(z1, z0, branch0, nu, num_scaled)
  num = val_branch$f
  #
  den_scaled = abs(x) > 400
  mybesselI = ifelse(abs(x) < 1000, Bessel::BesselI, Bessel::besselIasym)
  den = mybesselI(x, nu, expon.scaled = den_scaled)
  # if (!is.finite(den)) {
  #   fun = ifelse(abs(x) < 1000, "BesselI", "besselIasym")
  #   temp = "%s(%.4f, %.4f, expon.scaled=%d) = %.4f"
  #   stop(sprintf(temp, fun, x, nu, scaled, den))
  # }
  #
  f = f * (num/den)
  f = f * exp(num_scaled * abs(Re(z1)) -  den_scaled * abs(x))
  #
  return(list(f=f, z1=z1, branch1=val_branch$branch1))
}
