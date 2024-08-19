#' @rdname piv
#' @export
div <- function(x, v0, v1, tau, k, theta, sigma) {
  hj = seq(0.001, 1500, 0.001); j = 1:length(hj); h = hj[2] - hj[1]
  zs = CF(hj, v0, v1, tau, k, theta, sigma)
  #
  N = length(x)
  if (N == 1) {
    d = h/pi + (2/pi) * sum(cos(hj*x) * h * Re(zs))
  } else {
    d = rep(0, N)
    for (i in 1:N) {
      d[i] = h/pi + (2/pi) * sum(cos(hj*x[i]) * h * Re(zs))
    }
  }
  return(d)
}
