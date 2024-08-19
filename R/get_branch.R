#' Determine the correct branch
#'
#' @description
#' To accomplish the continuity of the Characteristic Function, we need to
#' track of arg(z) and change the branch when necessary, i.e., when the
#' complex numbers cross the negative real axis. see [CF()] and
#' [cBesselI()].
#' - `get_branch()` gets the branches all at once.
#' - `get_branch_next()` gets the branch one by one.
#' @name get_branch
#'
#' @param zs vector of complex numbers.
#' @param z1 current complex number.
#' @param z0 previous complex number.
#' @param branch0 branch of z0.
#'
#' @return
#' - `get_branch()`: vector of floors of the branches.
#' - `get_branch_next()`: scalar of floor of the current branch.
#' @export
#'
#' @examples
#' # a = seq(0, 1500, 0.001)
#' a = seq(0, 1500, 100)
#' k = 6.21; theta = 0.019; sigma = 0.61; v0 = 0.010201; v1 = v0; tau = 1
#' ga = sqrt(k^2 - 2*sigma^2*a*1i)
#' dk = exp(-k*tau); da = exp(-ga*tau)
#' zs = sqrt(v0*v1) * (4*ga/sigma^2) * exp(-0.5*ga*tau)/(1-da)
#' par(mfrow = c(1,1))
#' plot(Re(zs), Im(zs), col="blue")
#' branches = get_branch(zs)
#' #
#' z1 = -1 - 0.1i; z0 = -1 + 0.1i; branch0 = 0
#' branch1 = get_branch_next(z1, z0, branch0)
get_branch <- function(zs) {
  N = length(zs); branches = rep(0, N); floor = 0
  for (n in 1:N) {
    z  = zs[n]; if (n == 1) {next} # ignore the first one
    z0 = zs[n-1] # cross the negative real axis: increase the floor
    if (Arg(z) < 0 && Arg(z0) > 0) {floor = floor + 2*pi}
    branches[n] = floor            # arg = Arg(z) + floor
  }
  return(branches)
}
