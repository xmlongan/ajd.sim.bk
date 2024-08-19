#' @rdname cBesselI
#' @export
cBesselI_next <- function(z1, z0, branch0, v, expon_scaled=FALSE) {
  branch1 = get_branch_next(z1, z0, branch0)
  mybesselI = ifelse(abs(Re(z1)) < 1000, Bessel::BesselI, Bessel::besselIasym)
  cv = mybesselI(z1, v, expon.scaled = expon_scaled)
  # if (!is.finite(cv)) {
  #   fun = ifelse(abs(Re(z1)) < 1000, "BesselI", "besselIasym")
  #   temp = "%s(%.4f+%.4fi, %.4f, expon.scaled=%d) = %.4f+%.4fi"
  #   stop(sprintf(temp, fun, Re(z1), Im(z1), v, expon_scaled, Re(cv), Im(cv)))
  # }
  f = exp(v*branch1*1i) * cv
  return(list(f=f, branch1=branch1))
}
