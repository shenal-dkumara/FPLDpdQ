# Objective function: Hellinger distance between empirical and theoretical pdQ to find the shape parameters
.obj_shape <- function(par,u,fQstarhat,J){

  # if invalid parameters, return a large value
  if (!gld::gl.check.lambda(c(0, 1, par[1], par[2],par[3]), param = "fm5", vect = TRUE)) {
    return(.Machine$integer.max)
  }

  fQ <- gld::dqgl(u, c(0, 1, par[1], par[2],par[3]),param = "fm5")

  # Safely calculate normalising constant
  kappa <- tryCatch(
    integrate(gld::dqgl, 0, 1,
              lambda1 = 0, lambda2 = 1,
              lambda3 = par[1], lambda4 = par[2],
              lambda5 = par[3], param = "fm5")$value,
    error = function(e) NA_real_
  )

  if (!is.finite(kappa) || kappa <= 0) {
    return(.Machine$integer.max)
  }

  fQstar <- fQ / kappa
  integrand <- (sqrt(fQstarhat) - sqrt(fQstar)) ^ 2
  sqrt((sum(integrand) / J) / 2)
}
