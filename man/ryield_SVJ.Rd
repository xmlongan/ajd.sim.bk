% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ryield_SVJ.R
\name{ryield_SVJ}
\alias{ryield_SVJ}
\title{Simulate Yield using Broadie-Kaya Method for the SVJ model}
\usage{
ryield_SVJ(
  n,
  v0,
  tau,
  r,
  k,
  theta,
  sigma,
  rho,
  lambda,
  mu_bar,
  sigma_s,
  log.inv = F
)
}
\arguments{
\item{n}{number of yield samples to simulate.}

\item{v0}{current variance level.}

\item{tau}{time difference, \eqn{t-u}.}

\item{r}{riskless rate, risk-neutral drift.}

\item{k}{parameter \eqn{k}.}

\item{theta}{parameter \eqn{\theta}.}

\item{sigma}{parameter \eqn{\sigma}.}

\item{rho}{parameter \eqn{\rho}.}

\item{lambda}{parameter \eqn{\lambda}, arrival rate of jumps in the
underling asset price process.}

\item{mu_bar}{parameter \eqn{\bar{\mu}}.}

\item{sigma_s}{parameter \eqn{\sigma_s}.}

\item{log.inv}{\code{FALSE}, logging the Fourier inversion if \code{log.inv = TRUE}.}
}
\value{
next period yield, a scalar.
}
\description{
Simulate the next period yield given current variance using
Broadie-Kaya (2006) method for the SVJ model.
}
\examples{
v0 = 0.008836; k = 3.99; theta = 0.014; sigma = 0.27; rho = -0.79
r = 0.0319; tau = 5; lambda = 0.11; mu_bar = -0.12; sigma_s = 0.15
# Y = ryield_SVJ(1000, v0, tau, r, k, theta, sigma, rho,
# lambda, mu_bar, sigma_s)
# hist(Y)
}
