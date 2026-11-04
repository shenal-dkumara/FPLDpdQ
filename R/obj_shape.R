# Objective function: Hellinger distance between empirical and theoretical pdQ to find the FPLD shape parameters
.obj_shape <- function(par,u,fQstarhat,J){

  # Penalise invalid FPLD parameter combinations
  if (!gld::gl.check.lambda(c(0, 1, par[1], par[2],par[3]), param = "fm5", vect = TRUE)) {
    return(.Machine$integer.max)
  }

  # Theoretical quantile density
  fQ <- gld::dqgl(u, c(0, 1, par[1], par[2],par[3]),param = "fm5")

  # Calculate the normalising constant for the theoretical pdQ
  kappa <- tryCatch(
    integrate(gld::dqgl, 0, 1,
              lambda1 = 0, lambda2 = 1,
              lambda3 = par[1], lambda4 = par[2],
              lambda5 = par[3], param = "fm5")$value,
    error = function(e) NA_real_
  )

  # Penalise invalid normalising constants
  if (!is.finite(kappa) || kappa <= 0) {
    return(.Machine$integer.max)
  }
  # Normalise the theoretical quantile density to obtain the pdQ
  fQstar <- fQ / kappa

  # Approximate the Hellinger distance between empirical and theoretical pdQs
  integrand <- (sqrt(fQstarhat) - sqrt(fQstar)) ^ 2
  sqrt((sum(integrand) / J) / 2)
}
