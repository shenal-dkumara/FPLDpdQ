
#' Fit the Five-Parameter Lambda Distribution Using pdQ Methods
#'
#' @description
#' Estimates the parameters of the five-parameter lambda
#' distribution (FPLD) using probability density quantile
#' (pdQ) based methods.
#'
#' The shape parameters are estimated by minimising the
#' Hellinger distance between empirical and theoretical pdQ
#' functions. Location and scale parameters are subsequently
#' estimated using the selected method.
#'
#' @param x Numeric vector of observations.
#' @param J Number of grid points used for pdQ estimation.
#'   Default is 50.
#' @param method Method for estimating location and scale.
#'   One of "dens", "qm", "mom", or "qr".
#'   Default is "dens".
#' @param tau Quantile regression level, used when
#'   method = "qr". Default is 0.5.
#' @param weights Logical. If TRUE, weighted quantile
#'   regression is used when method = "qr".
#'   If FALSE, unweighted quantile regression is used.
#'   Default is TRUE.
#'
#' @return A named list containing the five estimated FPLD
#'   parameters: lambda1 (location), lambda2 (inverse scale),
#'   and lambda3, lambda4, and lambda5 (shape parameters).
#'
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#'
#' \donttest{
#' fit <- FPLDpdQ(x, method = "dens")
#' fit
#'
#' fit_wqr <- FPLDpdQ(x, method = "qr", weights = TRUE)
#' fit_wqr
#' }
#' @importFrom stats coef dcauchy integrate nlminb qcauchy quantile sd
#' @export

FPLDpdQ <- function(x, J=50, method="dens", tau=0.5, weights=TRUE){

  #remove NA values
  x <- x[!is.na(x)]

  # Empirical pdQ estimation
  inc <- 1/J
  u <- c(seq(inc / 2, 1 - inc / 2, inc))
  fQstarhat <- .pdQ_emp(x, J)

  # ---- shape levels (as starship) ----
  shape_levels <- list(
    lambda3 = c(-1.5, -1, -0.5, -0.1, 0, 0.1,0.2, 0.4, 0.8, 1, 1.5),
    lambda4 = c(-1.5, -1, -0.5, -0.1, 0, 0.1,0.2, 0.4, 0.8, 1, 1.5),
    lambda5 = c(-0.5, -0.25, 0, 0.25, 0.5)
  )

  # lower and upper bounds for shape parameters
  lower_shape <- c(-Inf, -Inf, -1)
  upper_shape <- c(Inf, Inf, 1)

  # lower and upper bounds for location and scale
  lower_ls <- c(-Inf, 0)
  upper_ls <- c(Inf, Inf)

  # ---- grid search for initial shape parameters ----
  grid_vals <- NMOF::gridSearch(fun = .obj_shape, levels = shape_levels, u=u, fQstarhat=fQstarhat, J=J,printDetail = FALSE)
  start_par <- grid_vals$minlevels

  # ---- optimize to find the optimal shape parameters ----
  opt_shape <- tryCatch(nlminb(start_par, .obj_shape, u = u, fQstarhat = fQstarhat, J = J,lower = lower_shape, upper = upper_shape ),error = function(e) list(par = c(NA,NA,NA)))
  shape_hat <- opt_shape$par
  l3 <- shape_hat[1]
  l4 <- shape_hat[2]
  l5 <- shape_hat[3]

  # ---- different methods for estimating location and scale parematers ----
  if(method=="qm"){

    # quartile matching estimates for l1 and l2
    pr<-c(0.25,0.5,0.75)
    q<-c(quantile(x,pr[1]),quantile(x,pr[2]),quantile(x,pr[3]))

    # Handle zero shape parameters
    term1 <- if (l3 == 0) log(pr) else expm1(l3 * log(pr)) / l3
    term2 <- if (l4 == 0) log(1 - pr) else expm1(l4 * log(1 - pr)) / l4

    xpr <- (1 - l5) * term1 - (1 + l5) * term2

    l2 <- (xpr[3]-xpr[1])/(q[3]-q[1])
    l1 <- q[2]-xpr[2]/l2
    est_par <- c(l1,l2,l3,l4,l5)

  }else if(method=="dens"){

    ## ---- Start values from QM (unchanged) ----
    pr <- c(0.25,0.5,0.75)
    q  <- quantile(x, pr, names = FALSE)

    # Handle zero shape parameters
    term1 <- if (l3 == 0) log(pr) else expm1(l3 * log(pr)) / l3
    term2 <- if (l4 == 0) log(1 - pr) else expm1(l4 * log(1 - pr)) / l4
    xpr <- (1 - l5) * term1 - (1 + l5) * term2

    l2.start <- (xpr[3]-xpr[1])/(q[3]-q[1])
    l1.start <- q[2] - xpr[2]/l2.start

    start_ls <- c(l1.start, l2.start)

    ## ---- PRECOMPUTE EMPIRICAL DENSITY ONCE ----
    n  <- length(x)
    emp <- EnvStats::demp(x, x)
    sqrt_emp <- sqrt(pmax(emp, 0))

    ## ---- Fast objective ----
    obj_dens <- function(par){

      # Fast reject
      if (!gld::gl.check.lambda(c(par[1], par[2], l3, l4, l5),
                                param = "fm5", vect = TRUE)) {
        return(.Machine$integer.max)
      }

      th <- gld::dgl(x, par[1], par[2], l3, l4, l5, param = "fm5")
      sqrt_th <- sqrt(pmax(th, 0))

      # direct mean is faster than sum/length
      sqrt(0.5 * mean((sqrt_emp - sqrt_th)^2))
    }

    opt_ls <- tryCatch(
      nlminb(start_ls, obj_dens,
             lower = lower_ls,
             upper = upper_ls),
      error = function(e) list(par = c(NA,NA))
    )

    l1 <- opt_ls$par[1]
    l2 <- opt_ls$par[2]
    est_par <- c(l1,l2,l3,l4,l5)
  }else if (method == "mom") {

    l1 <- NA_real_
    l2 <- NA_real_

    p_grid <- seq(1e-4, 1 - 1e-4, length.out = 2000)

    try({S_vals <- gld::qgl(p_grid,lambda1 = c(0, 1, l3, l4, l5),param = "fm5")

    S_vals <- S_vals[is.finite(S_vals)]

    if (length(S_vals) > 10) {

      S1 <- mean(S_vals)
      S2 <- mean(S_vals^2)

      var_S <- S2 - S1^2

      if (is.finite(S1) && is.finite(var_S) && var_S > 0 && sd(x) > 0) {

        # lambda2 is inverse scale
        l2 <- sqrt(var_S) / sd(x)

        # lambda1 = sample mean - S1 / lambda2
        l1 <- mean(x) - S1 / l2
      }
    }

    }, silent = TRUE)

    est_par <- c(l1, l2, l3, l4, l5)

} else if (method == "qr") {

  n <- length(x)
  pr <- (1:n) / (n + 1)
  x_ord <- sort(x)

  # Theoretical quantiles using estimated shape parameters
  qgl_theo <- gld::qgl(
    pr,
    c(0, 1, l3, l4, l5),
    param = "fm5"
  )

  if (weights) {

    # Weighted quantile regression (pdQ-WQR)

    # Theoretical density evaluated at quantiles
    fQ <- gld::dqgl(
      pr,
      c(0, 1, l3, l4, l5),
      param = "fm5"
    )

    # Inverse asymptotic quantile SE weights
    w <- fQ / sqrt(pr * (1 - pr))

    # Numerical safeguards
    w[!is.finite(w) | w < 0] <- 0

    if (!any(w > 0)) {
      stop("No positive finite QR weights")
    }

    # Normalise weights
    w <- w / max(w)

    # Weighted quantile regression
    qr_fit <- quantreg::rq(
      x_ord ~ qgl_theo,
      tau = tau,
      weights = w,
      method = "fn"
    )

    l1 <- unname(coef(qr_fit)[1])
    slope <- unname(coef(qr_fit)[2])

  } else {

    # Unweighted quantile regression (pdQ-QR)
    X <- cbind(1, qgl_theo)

    qr_fit <- quantreg::rq.fit(
      x = X,
      y = x_ord,
      tau = tau,
      method = "fn"
    )

    l1 <- qr_fit$coefficients[1]
    slope <- qr_fit$coefficients[2]
  }

  if (!is.finite(slope) || slope <= 0) {
    stop("Invalid QR slope")
  }

  l2 <- 1 / slope

  est_par <- c(l1, l2, l3, l4, l5)

} else{
  stop("method must be one of 'qm', 'dens', 'mom', or 'qr'")
}
  names(est_par) <- c("lambda1","lambda2","lambda3","lambda4","lambda5")
  as.list(est_par)
}
