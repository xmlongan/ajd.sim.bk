#' @rdname get_branch
#' @export
get_branch_next <- function(z1, z0, branch0) {
  # if (any(!is.finite(c(z1, z0)))) {
  #   temp = "z1=(%.7f, %.7f), z0=(%.7f, %.7f)"
  #   stop(sprintf(temp, Re(z1), Im(z1), Re(z0), Im(z0)))
  # }
  if (Arg(z1) < 0 && Arg(z0) > 0) { return(branch0 + 2*pi) }
  else {return(branch0)}
}
