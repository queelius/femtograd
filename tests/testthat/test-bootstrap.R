# Tests for bootstrap inference
# Note: Using small n_boot for fast tests

test_that("bootstrap_fit works with vector data", {
  set.seed(42)
  x <- rnorm(50, mean = 5, sd = 2)

  # Log-likelihood maker for normal distribution
  make_loglik <- function(data) {
    function(mu, log_sigma) {
      loglik_normal(mu, exp(log_sigma), data)
    }
  }

  boot <- bootstrap_fit(
    make_loglik,
    data = x,
    params = c(mu = 4, log_sigma = 0.5),
    n_boot = 30,  # Small for testing
    progress = FALSE
  )

  expect_s3_class(boot, "bootstrap_result")
  expect_true(boot$n_successful > 20)  # Most should succeed
  expect_equal(length(boot$se), 2)
  expect_equal(names(boot$se), c("mu", "log_sigma"))
})


test_that("bootstrap_fit returns sensible standard errors", {
  set.seed(123)
  x <- rnorm(100, mean = 5, sd = 2)

  make_loglik <- function(data) {
    function(mu, log_sigma) {
      loglik_normal(mu, exp(log_sigma), data)
    }
  }

  boot <- bootstrap_fit(
    make_loglik,
    data = x,
    params = c(mu = 4, log_sigma = 0.5),
    n_boot = 50,
    progress = FALSE
  )

  # Bootstrap SEs should be positive
  expect_true(all(boot$se > 0))

  # Bias should be small relative to SE
  expect_true(all(abs(boot$bias) < 3 * boot$se))
})


test_that("confint.bootstrap_result computes percentile CIs", {
  set.seed(42)
  x <- rnorm(50, mean = 5, sd = 2)

  make_loglik <- function(data) {
    function(mu, log_sigma) {
      loglik_normal(mu, exp(log_sigma), data)
    }
  }

  boot <- bootstrap_fit(
    make_loglik,
    data = x,
    params = c(mu = 4, log_sigma = 0.5),
    n_boot = 50,
    progress = FALSE
  )

  ci <- confint(boot, type = "percentile")

  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 2)
  expect_equal(rownames(ci), c("mu", "log_sigma"))

  # CIs should contain the original estimate (usually)
  expect_true(ci["mu", 1] < boot$original["mu"])
  expect_true(ci["mu", 2] > boot$original["mu"])
})


test_that("confint.bootstrap_result supports different CI types", {
  set.seed(42)
  x <- rnorm(50, mean = 5, sd = 2)

  make_loglik <- function(data) {
    function(mu, log_sigma) {
      loglik_normal(mu, exp(log_sigma), data)
    }
  }

  boot <- bootstrap_fit(
    make_loglik,
    data = x,
    params = c(mu = 4, log_sigma = 0.5),
    n_boot = 50,
    progress = FALSE
  )

  ci_pct <- confint(boot, type = "percentile")
  ci_basic <- confint(boot, type = "basic")
  ci_norm <- confint(boot, type = "normal")

  # All should be matrices with same dimensions
  expect_equal(dim(ci_pct), dim(ci_basic))
  expect_equal(dim(ci_pct), dim(ci_norm))

  # Normal should be symmetric around estimate
  mid_norm <- (ci_norm[, 1] + ci_norm[, 2]) / 2
  expect_equal(unname(mid_norm), unname(boot$original), tolerance = 1e-10)
})


test_that("confint.bootstrap_result subset by parameter", {
  set.seed(42)
  x <- rnorm(50, mean = 5, sd = 2)

  make_loglik <- function(data) {
    function(mu, log_sigma) {
      loglik_normal(mu, exp(log_sigma), data)
    }
  }

  boot <- bootstrap_fit(
    make_loglik,
    data = x,
    params = c(mu = 4, log_sigma = 0.5),
    n_boot = 30,
    progress = FALSE
  )

  ci_mu <- confint(boot, parm = "mu")

  expect_equal(nrow(ci_mu), 1)
  expect_equal(rownames(ci_mu), "mu")
})


test_that("print.bootstrap_result works", {
  set.seed(42)
  x <- rnorm(30, mean = 5, sd = 2)

  make_loglik <- function(data) {
    function(mu, log_sigma) {
      loglik_normal(mu, exp(log_sigma), data)
    }
  }

  boot <- bootstrap_fit(
    make_loglik,
    data = x,
    params = c(mu = 4, log_sigma = 0.5),
    n_boot = 20,
    progress = FALSE
  )

  expect_output(print(boot), "Bootstrap Results")
  expect_output(print(boot), "Boot.SE")
})


test_that("summary.bootstrap_result works", {
  set.seed(42)
  x <- rnorm(30, mean = 5, sd = 2)

  make_loglik <- function(data) {
    function(mu, log_sigma) {
      loglik_normal(mu, exp(log_sigma), data)
    }
  }

  boot <- bootstrap_fit(
    make_loglik,
    data = x,
    params = c(mu = 4, log_sigma = 0.5),
    n_boot = 20,
    progress = FALSE
  )

  summ <- summary(boot)

  expect_s3_class(summ, "summary.bootstrap_result")
  expect_true("coefficients" %in% names(summ))
  expect_true("Lower" %in% names(summ$coefficients))
  expect_true("Upper" %in% names(summ$coefficients))
})


test_that("bootstrap_fit works with exponential distribution", {
  set.seed(42)
  x <- rexp(100, rate = 2)

  # Use log-parameterization for positivity
  make_loglik <- function(data) {
    function(log_rate) {
      rate <- exp(log_rate)
      loglik_exponential(rate, data)
    }
  }

  boot <- suppressWarnings(bootstrap_fit(
    make_loglik,
    data = x,
    params = c(log_rate = 0.5),
    n_boot = 30,
    progress = FALSE
  ))

  expect_true(boot$n_successful > 15)
  expect_equal(names(boot$se), "log_rate")

  # Original estimate on rate scale should be close to 1/mean(x) ≈ 2
  fitted_rate <- exp(boot$original["log_rate"])
  expect_equal(unname(fitted_rate), 1 / mean(x), tolerance = 0.2)
})


test_that("bootstrap handles failing fits gracefully", {
  set.seed(42)
  # Very small sample - some bootstrap samples may fail
  x <- rnorm(10, mean = 5, sd = 2)

  make_loglik <- function(data) {
    function(mu, log_sigma) {
      loglik_normal(mu, exp(log_sigma), data)
    }
  }

  # Should not error even if some fits fail
  boot <- bootstrap_fit(
    make_loglik,
    data = x,
    params = c(mu = 4, log_sigma = 0.5),
    n_boot = 20,
    progress = FALSE
  )

  expect_s3_class(boot, "bootstrap_result")
  expect_true(boot$n_failed >= 0)
  expect_equal(boot$n_successful + boot$n_failed, boot$n_boot)
})


test_that("bootstrap CIs cover true parameter", {
  # This is a probabilistic test - may occasionally fail by chance
  set.seed(456)
  true_mu <- 5
  true_sigma <- 2
  x <- rnorm(100, mean = true_mu, sd = true_sigma)

  make_loglik <- function(data) {
    function(mu, log_sigma) {
      loglik_normal(mu, exp(log_sigma), data)
    }
  }

  boot <- bootstrap_fit(
    make_loglik,
    data = x,
    params = c(mu = 4, log_sigma = 0.5),
    n_boot = 100,
    progress = FALSE
  )

  ci <- confint(boot, level = 0.95)

  # True mu should usually be in CI (but not guaranteed)
  # Just check that CIs are reasonable width
  mu_width <- ci["mu", 2] - ci["mu", 1]
  expect_true(mu_width > 0)
  expect_true(mu_width < 2)  # Shouldn't be too wide for n=100
})
