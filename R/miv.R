#' First Two Conditional Moments of IV
#'
#' @description
#' Compute the first two conditional moments of IV, given \eqn{v_u} and
#'  \eqn{v_t}.
#'
#' @param v0 the left side variance level, \eqn{v_u}.
#' @param v1 the right side variance level, \eqn{v_t}.
#' @param tau time difference, \eqn{t-u}.
#' @param k parameter \eqn{k}.
#' @param theta parameter \eqn{\theta}.
#' @param sigma parameter \eqn{\sigma}.
#'
#' @return vector of the first two un-centralized moments \eqn{\mu_1, \mu_2}.
#' @export
#'
#' @examples
#' k = 6.21; theta = 0.019; sigma = 0.61; v0 = 0.010201; v1 = v0; tau = 1
#' mu = miv(v0, v1, tau, k, theta, sigma)
#' sd = sqrt(mu[2] - mu[1]^2)
#' x = seq(mu[1]-4*sd, mu[1]+4*sd, 0.001)
#' plot(x, dnorm(x, mu[1], sd), type="l", col="blue")
miv <- function(v0, v1, tau, k, theta, sigma) {
  d = 4 * theta * k / sigma^2
  v = 0.5*d - 1
  #
  dk = exp(-k*tau)
  x = sqrt(v0*v1) * 4*k*exp(-0.5*k*tau) / (sigma^2*(1-dk))
  #
  scaled = abs(x) > 400
  mybesselI = ifelse(abs(x) < 1000, Bessel::BesselI, Bessel::besselIasym)
  Ix  = mybesselI(x, v,   expon.scaled = scaled)
  nv1 = mybesselI(x, v-1, expon.scaled = scaled)
  nv2 = mybesselI(x, v+1, expon.scaled = scaled)
  nv3 = mybesselI(x, v-2, expon.scaled = scaled)
  nv4 = mybesselI(x, v+2, expon.scaled = scaled)
  # if (any(!is.finite(c(Ix, nv1, nv2, nv3, nv4)))) {
  #   fun = ifelse(abs(x) < 1000, "BesselI", "besselIasym")
  #   temp = "%s(%.4f, %.4f, expon.scaled=%d) = %.4f"
  #   stop(sprintf(temp, fun, x, v, scaled, Ix))
  # }
  #
  dI_I = nv1/(2*Ix) + nv2/(2*Ix)
  ddI_I= nv3/(4*Ix) + Ix/(2*Ix) + nv4/(4*Ix)
  #
  df1 = (sigma^2/k)*(-1/k + tau/2 + tau*dk/(1-dk)) # times i
  df2 = (v0+v1) * ((1/k)*(1+dk)/(1-dk) - 2*tau*dk/(1-dk)^2) # times i
  df3 = (-4/k+2*tau)*exp(-0.5*k*tau) + 4*tau*exp(-1.5*k*tau)/(1-dk)
  df3 = (df3/(1-dk)) * sqrt(v0*v1) * dI_I # times i
  #
  mu1 = df1 + df2 + df3
  #
  ddf1 = 1/k^2 + tau/(2*k) - tau^2/4 + (tau/k - 2*tau^2)*dk/(1-dk)
  ddf1 = ddf1 - 2*tau^2*dk^2/(1-dk)^2
  ddf1 = ddf1 * (sigma^4/k^2)
  #
  t1 = (1+dk)/k^2 + tau*dk/k - tau^2*dk
  t2 = (tau*(1+dk)*dk/k - 3*tau^2*dk^2 - tau^2*dk) / (1-dk)
  t3 = 2*tau^2*dk^2*(1+dk)/(1-dk)^2
  ddf2 = - df2^2 - ((v0+v1)*sigma^2/(k*(1-dk))) * (t1 + t2 - t3)
  #
  dz = (-4/k+2*tau)*exp(-0.5*k*tau) + 4*tau*exp(-1.5*k*tau)/(1-dk)
  dz = dz * sqrt(v0*v1)/(1-dk) # times i
  #
  t1 = (4/k^2 + 2*tau/k - tau^2) * exp(-0.5*k*tau)
  t2 = (4*tau/k - 8*tau^2) * exp(-1.5*k*tau) / (1-dk)
  t3 = 8*tau^2 * exp(-2.5*k*tau) / (1-dk)^2
  ddz = sqrt(v0*v1)*sigma^2/(k*(1-dk)) * (t1 + t2 - t3)
  ddf3 = - dz^2 * ddI_I + dI_I * ddz
  #
  mu2 = ddf1 + ddf2 + ddf3 - 2*(df1*df2 + df1*df3 + df2*df3)
  mu2 = - mu2
  #
  return(c(mu1, mu2))
}
