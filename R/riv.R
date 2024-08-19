#' @rdname piv
riv <- function(n, v0, v1, tau, k, theta, sigma, log.inv = FALSE) {
  #                       Cautions
  # when tau is close to 0, it will produce large x or z
  # for the modified Bessel function of the first kind I_v(x).
  # Consequently, approxiamtion errors would probabily occur
  # if we keep the normal procedure.
  #
  #                   Special Preprocessing
  # when tau is very small, assume IV is fixed
  if (tau < 0.0001) {
    # temp = "time_step = %.7f, too small, set iv = tau*(v0+v1)/2\n"
    # cat(sprintf(temp, tau));
    return(rep(tau*(v0+v1)/2, n))
  }
  #                Set up the grid size h
  # through the first two moments
  mu = miv(v0, v1, tau, k, theta, sigma)
  #
  if (any(!is.finite(mu))) stop(sprintf("mu[1]=%.20f, mu[2]=%.20f",mu[1],mu[2]))
  var = mu[2] - (mu[1])^2
  if (var <= 0) {
    # temp = "miv(%.10f, %.10f, %.10f, %.4f, %.4f, %.4f)"
    # msg1 = sprintf(temp, v0, v1, tau, k, theta, sigma)
    # temp = "var(iv) = %.20f, mu[1]=%.20f, mu[2]=%.20f"
    # msg2 = sprintf(temp, var, mu[1], mu[2])
    # # stop(paste("var(iv) <= 0", msg1, msg2, sep="\n"))
    # cat(paste("var(iv) <= 0", msg1, msg2, "", sep="\n"))
    return(rep(tau*(v0+v1)/2, n))
  }
  #
  sd = sqrt(var)
  mu_eps = mu[1] + 7*sd; h = pi/mu_eps
  #
  #     Truncate the Characteristic Function inversion
  # through error bound
  #
  if (log.inv) { t.start = Sys.time() }
  #
  zs = CF_trunc(h, v0, v1, tau, k, theta, sigma, eps = 1e-5)
  j    = seq(1, length(zs))
  hj   = h*j
  Rezs = Re(zs); rm(zs)
  #
  if (log.inv) {
    t.taken = Sys.time() - t.start; t.taken = as.numeric(t.taken)
    fname = format(t.start, "./script/logs/riv-%Hhour-%Mmin.csv")
    if (!file.exists(fname)) {
      title = c("v0", "v1", "tau", "k", "theta", "sigma", "h",
                "num_Bessel_eval", "secs_consumed\n")
      cat(title, file=fname, sep=',', append=T)
    }
    str = paste(v0, v1, tau, k, theta, sigma, h, 2*length(j), t.taken, sep=',')
    cat(str, '\n', file=fname, append=T)
  }
  #
  xs = rep(0, n)
  for (k in 1:n) {
    U = stats::runif(1)
    #             normal approximation
    x0 = stats::qnorm(U, mean=mu[1], sd=sd)
    x0 = ifelse(x0 > 0, x0, 0.01 * mu[1])       # keep positive
    #          second-order Newton's method
    # sovle f(x) = piv(x, k, theta, sigma, v0, v1, tau) - U = 0
    iter = 0
    repeat {
      iter = iter + 1
      f   = h*x0/pi + (2/pi) * sum((sin(hj*x0)/j) * Rezs) - U
      df  = h/pi    + (2/pi) * sum(cos(hj*x0) * h * Rezs)
      ddf =         - (2/pi) * sum(sin(hj*x0) * hj * h * Rezs)
      delta = 1 - 2*f*ddf/df^2
      if (delta >= 0) { x = x0 - (df/ddf) * (1 - sqrt(delta)) }
      else { found = FALSE; break }
      if (abs(x-x0) < 1e-5 || iter > 10) { found = TRUE; break }
      x0 = x
    }
    if (!found || x <= 0) {
      # cat(sprintf("%d-iter: resort to bisect, delta = %.7f\n", iter, delta))
      #           resort to bisection
      lb = 0; ub = mu[1] + 7*sd; x = (lb + ub)/2
      while (ub - lb > 1e-5) {
        x = (lb + ub)/2
        p = h*x/pi + (2/pi) * sum((sin(hj*x)/j) * Rezs)
        if (p < U)      { lb = x }
        else if (p > U) { ub = x }
        else            { break  } # found, break the while loop
      }
    }
    xs[k] = x
  }
  return(xs)
}
