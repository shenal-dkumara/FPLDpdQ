test_that("FPLDpdQ returns valid estimates for all methods", {

  set.seed(123)
  x <- rnorm(100)

  fits <- list(
    dens = FPLDpdQ(x, method = "dens"),
    qm   = FPLDpdQ(x, method = "qm"),
    mom  = FPLDpdQ(x, method = "mom"),
    qr   = FPLDpdQ(x, method = "qr", weights = FALSE),
    wqr  = FPLDpdQ(x, method = "qr", weights = TRUE)
  )

  for (fit in fits) {

    # Check that five parameters are returned
    expect_length(fit, 5)

    # Check parameter names
    expect_named(
      fit,
      c("lambda1", "lambda2", "lambda3",
        "lambda4", "lambda5")
    )

    # Check that all estimates are finite
    expect_true(all(is.finite(unlist(fit))))

    # Check that inverse scale is positive
    expect_gt(fit$lambda2, 0)
  }

  # Shape estimates should be identical across methods
  shapes <- lapply(fits, function(fit) {
    unlist(fit[c("lambda3", "lambda4", "lambda5")])
  })

  for (i in 2:length(shapes)) {
    expect_equal(shapes[[i]], shapes[[1]])
  }

})
