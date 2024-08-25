#' @rdname piv
riv <- function(v0, v1, tau, k, theta, sigma, log.inv = FALSE) {
  #                       Cautions
  # when tau is close to 0, it will produce large x or z
  # for the modified Bessel function of the first kind I_v(x).
  # Consequently, approximation errors would probably occur
  # if we keep the normal procedure.
  #                   Special Preprocessing
  # when tau is very small, assume IV is fixed
  if (tau < 0.0001) {
    # temp = "time_step = %.7f, too small, set iv = tau*(v0+v1)/2\n"
    # cat(sprintf(temp, tau));
    return(tau*(v0+v1)/2)
  }
  #               Compute the first two moments
  mu = miv(v0, v1, tau, k, theta, sigma)
  if (any(!is.finite(mu))) {
    stop(sprintf("mu[1]=%.20f, mu[2]=%.20f",mu[1],mu[2]))
  }
  var = mu[2] - (mu[1])^2
  if (var <= 0) {
    # temp = "miv(%.10f, %.10f, %.10f, %.4f, %.4f, %.4f)"
    # msg1 = sprintf(temp, v0, v1, tau, k, theta, sigma)
    # temp = "var(iv) = %.20f, mu[1]=%.20f, mu[2]=%.20f"
    # msg2 = sprintf(temp, var, mu[1], mu[2])
    # cat(paste("var(iv) <= 0", msg1, msg2, "", sep="\n"))
    return(tau*(v0+v1)/2)
  }
  # discretization error setting
  sd = sqrt(var); mu_eps = mu[1] + 5*sd
  #
  U = stats::runif(1)
  #                    Normal approximation
  x0 = stats::qnorm(U, mean=mu[1], sd=sd); x0_original = x0
  x0 = ifelse(x0 > 0, x0, 0.01 * mu[1])       # keep positive
  iter = 0
  repeat {
    iter = iter + 1
    #           Set up grid size to bound discretization error
    h = 2*pi/(x0 + mu_eps)
    #      Truncate the Characteristic Function inversion
    zs = CF_trunc(h, v0, v1, tau, k, theta, sigma, eps=1e-5, log.inv=log.inv)
    j = seq(1, length(zs)); hj = h*j; Rezs = Re(zs); rm(zs)
    #
    #              Second-order Newton's method
    # sovle f(x) = piv(x, k, theta, sigma, v0, v1, tau) - U = 0
    #
    f   = h*x0/pi + (2/pi) * sum((sin(hj*x0)/j) * Rezs) - U
    df  = h/pi    + (2/pi) * sum(cos(hj*x0) * h * Rezs)
    ddf =         - (2/pi) * sum(sin(hj*x0) * hj * h * Rezs)
    delta = 1 - 2*f*ddf/df^2
    if (delta >= 0) {
      x = x0 - (df/ddf) * (1 - sqrt(delta))
      if (x < 0) {found = FALSE; break }
    } else { found = FALSE; x = x0; break }
    if (abs(x-x0) < 1e-5 || iter > 10) { found = TRUE; break }
    x0 = x
  }
  if (!found) {#        Resort to bisection
    # cat(sprintf("iterated: %d times\n", iter))
    # cat(sprintf("  resort to bisect, delta = %.7f, x = %.7f\n", delta, x))
    LB = 0; UB = mu[1]+5*sd
    # temp = "  (LB, UB)=(%.7f, %.7f), mean=%.7f, sd=%.7f, U=%.7f, x0=%.7f\n"
    # cat(sprintf(temp, LB, UB, mu[1], sd, U, x0_original))
    #
    # Set up lower bound and upper bound
    #
    if (U < 0.1) {
      lb = LB; ub = mu[1] - sd
      if (ub < 0) { n_sd = mu[1]/sd; ub = mu[1] - (n_sd/2) * sd }
    }
    else if (U > 0.9) { lb = mu[1] + 2*sd; ub = UB }
    else { lb = max(x0_original-2*sd, LB); ub = min(x0_original+2*sd, UB) }
    # cat(sprintf("  (lb, ub) = (%.7f, %.7f)\n", lb, ub))
    #
    x = (lb + ub)/2 # lb = 0; ub = mu[1] + 5*sd; x = (lb + ub)/2
    # i = 0
    while (ub - lb > 1e-5) {
      # i = i + 1
      x = (lb + ub)/2
      # cat(sprintf("    (lb, ub) = (%.7f, %.7f), x = %.7f\n", lb, ub, x))
      h = 2*pi/(x + mu_eps)
      zs = CF_trunc(h, v0, v1, tau, k, theta, sigma, eps=1e-5, log.inv=log.inv)
      j = seq(1, length(zs)); hj = h*j; Rezs = Re(zs); rm(zs)
      #
      p = h*x/pi + (2/pi) * sum((sin(hj*x)/j) * Rezs)
      if (abs(p-U) < 1e-5) { break }  # found, break the while loop
      else if (p < U)      { lb = x }
      else                 { ub = x }
    }
    # cat(sprintf("  bisected: %d times\n", i))
  }
  return(x)
}
