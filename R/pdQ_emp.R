# Estimate the quantile density using an Epanechnikov kernel
.quantile_dens <- function(u, qratio, x) {
  n <- length(x)

  # Calculate the bandwidth
  Ap <- (15 ^ (1 / 5)) * abs(qratio) ^ (2 / 5)
  bw <- Ap / (n ^ (1 / 5))

  # Sort observations
  xsort <- sort(x)

  # Epanechnikov kernel function
  kernepach <- function(u) (1 - u ^ 2) * (abs(u) <= 1) * 3 / 4

  # Kernel-based quantile density estimate
  consts <- kernepach((u - (1:n - 1)/n)/bw)/bw -
    kernepach((u - (1:n)/n)/bw)/bw
  sum(xsort * consts)
}

# Estimate the empirical probability density quantile (pdQ) function
.pdQ_emp <- function(x, J) {

  # Equally spaced midpoint grid
  inc <- 1 / J
  u <- seq(inc / 2, 1 - inc / 2, by = inc)

  # Calculate bandwidth components using the Cauchy reference distribution
  qcauchyu    <- 1 / dcauchy(qcauchy(u))
  theta       <- pi * (u - 0.5)
  qcauchypr   <- 2 * pi^2 * (tan(theta) + tan(theta)^3)
  qcauchyprpr <- 2 * pi^3 * (1 + 3 * tan(theta)^2) / (cos(theta)^2)
  qratio <- qcauchyu / qcauchyprpr

  # Estimate the quantile density
  qhat <- mapply(.quantile_dens, u, qratio, MoreArgs = list(x = x))

  # Obtain the empirical density quantile by taking the reciprocal
  fQhat <- 1 / qhat

  # Approximate the normalising constant
  kappahat <- inc * sum(fQhat)

  # Normalise to obtain the empirical pdQ
  fQstarhat <- fQhat / kappahat
  fQstarhat
}
