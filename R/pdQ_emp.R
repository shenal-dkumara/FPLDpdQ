# Quantile Density Estimation
.quantile_dens <- function(u, qratio, x) {
  n <- length(x)
  Ap <- (15 ^ (1 / 5)) * abs(qratio) ^ (2 / 5)
  bw <- Ap / (n ^ (1 / 5))
  xsort <- sort(x)
  kernepach <- function(u) (1 - u ^ 2) * (abs(u) <= 1) * 3 / 4
  consts <- kernepach((u - (1:n - 1)/n)/bw)/bw -
    kernepach((u - (1:n)/n)/bw)/bw
  sum(xsort * consts)
}

# Empirical estimation of the pdQ
.pdQ_emp <- function(x, J) {
  inc <- 1 / J
  u <- seq(inc / 2, 1 - inc / 2, by = inc)
  qcauchyu    <- 1 / dcauchy(qcauchy(u))
  theta       <- pi * (u - 0.5)
  qcauchypr   <- 2 * pi^2 * (tan(theta) + tan(theta)^3)
  qcauchyprpr <- 2 * pi^3 * (1 + 3 * tan(theta)^2) / (cos(theta)^2)
  qratio <- qcauchyu / qcauchyprpr

  ## Estimate quantile density
  qhat <- mapply(.quantile_dens, u, qratio, MoreArgs = list(x = x))

  ## Compute fQ and normalize
  fQhat <- 1 / qhat
  kappahat <- inc * sum(fQhat)
  fQstarhat <- fQhat / kappahat
  fQstarhat
}
