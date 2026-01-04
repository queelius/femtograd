# Tests for fit() function and femtofit class

# Helper: Normal log-likelihood with log(sigma) parameterization
# This is the standard approach for constrained positive parameters
loglik_normal_logsig <- function(mu, log_sigma, x) {
  sigma <- exp(log_sigma)
  loglik_normal(mu, sigma, x)
}

test_that("fit() works with named arguments style", {
  set.seed(42)
  x <- rnorm(100, mean = 5, sd = 2)

  # Use log(sigma) parameterization to keep sigma > 0
  result <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)  # log(1) = 0
  )

  expect_s3_class(result, "femtofit")
  expect_true(result$converged)

  # Check estimates are close to true values
  expect_equal(unname(coef(result)["mu"]), mean(x), tolerance = 0.2)
  expect_equal(unname(exp(coef(result)["log_sigma"])), sd(x), tolerance = 0.3)
})


test_that("fit() works with single parameter argument style", {
  set.seed(42)
  x <- rnorm(100, mean = 5, sd = 2)

  result <- fit(
    function(p) loglik_normal_logsig(p$mu, p$log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  expect_s3_class(result, "femtofit")
  expect_true(result$converged)

  # Check estimates
  expect_equal(unname(coef(result)["mu"]), mean(x), tolerance = 0.2)
})


test_that("fit() works with different optimization methods", {
  set.seed(123)
  x <- rexp(50, rate = 2)

  # BFGS (default)
  result_bfgs <- fit(
    function(rate) loglik_exponential(rate, x),
    params = c(rate = 1),
    method = "bfgs"
  )
  expect_true(result_bfgs$converged)

  # L-BFGS
  result_lbfgs <- fit(
    function(rate) loglik_exponential(rate, x),
    params = c(rate = 1),
    method = "lbfgs"
  )
  expect_true(result_lbfgs$converged)

  # Newton-Raphson
  result_newton <- fit(
    function(rate) loglik_exponential(rate, x),
    params = c(rate = 1),
    method = "newton"
  )
  expect_true(result_newton$converged)

  # All should give similar estimates
  expect_equal(coef(result_bfgs), coef(result_lbfgs), tolerance = 0.01)
  expect_equal(coef(result_bfgs), coef(result_newton), tolerance = 0.01)
})


test_that("coef() returns named vector", {
  set.seed(42)
  x <- rnorm(50, mean = 3, sd = 1)

  result <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  cf <- coef(result)
  expect_type(cf, "double")
  expect_named(cf, c("mu", "log_sigma"))
})


test_that("vcov() returns named matrix", {
  set.seed(42)
  x <- rnorm(50, mean = 3, sd = 1)

  result <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  V <- vcov(result)
  expect_true(is.matrix(V))
  expect_equal(dim(V), c(2, 2))
  expect_equal(rownames(V), c("mu", "log_sigma"))
  expect_equal(colnames(V), c("mu", "log_sigma"))

  # Should be symmetric
  expect_equal(V, t(V))

  # Should be positive definite (for a well-identified model)
  expect_true(all(eigen(V)$values > 0))
})


test_that("confint() returns confidence intervals", {
  set.seed(42)
  x <- rnorm(100, mean = 5, sd = 2)

  result <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  # Default 95% CI
  ci <- confint(result)
  expect_true(is.matrix(ci))
  expect_equal(dim(ci), c(2, 2))
  expect_equal(rownames(ci), c("mu", "log_sigma"))

  # Lower < estimate < upper
  expect_true(all(ci[, 1] < coef(result)))
  expect_true(all(coef(result) < ci[, 2]))

  # 90% CI should be narrower
  ci90 <- confint(result, level = 0.90)
  expect_true(all((ci90[, 2] - ci90[, 1]) < (ci[, 2] - ci[, 1])))

  # Subset by parameter name
  ci_mu <- confint(result, parm = "mu")
  expect_equal(nrow(ci_mu), 1)
  expect_equal(rownames(ci_mu), "mu")
})


test_that("logLik() returns proper logLik object", {
  set.seed(42)
  x <- rnorm(50, mean = 3, sd = 1)

  result <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  ll <- logLik(result)
  expect_s3_class(ll, "logLik")
  expect_equal(attr(ll, "df"), 2)  # 2 parameters
})


test_that("AIC() and BIC() work via logLik()", {
  set.seed(42)
  x <- rnorm(50, mean = 3, sd = 1)

  result <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  # AIC should work
  aic_val <- AIC(result)
  expect_type(aic_val, "double")
  expect_true(is.finite(aic_val))

  # AIC = -2 * loglik + 2 * k
  expected_aic <- -2 * result$loglik + 2 * 2
  expect_equal(aic_val, expected_aic)
})


test_that("se() returns standard errors", {
  set.seed(42)
  x <- rnorm(100, mean = 5, sd = 2)

  result <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  ses <- se(result)
  expect_type(ses, "double")
  expect_named(ses, c("mu", "log_sigma"))
  expect_true(all(ses > 0))

  # SE should equal sqrt(diag(vcov))
  expect_equal(ses, sqrt(diag(vcov(result))))
})


test_that("summary() produces coefficient table", {
  set.seed(42)
  x <- rnorm(50, mean = 3, sd = 1)

  result <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  s <- summary(result)
  expect_s3_class(s, "summary.femtofit")

  # Coefficient table
  expect_true(is.matrix(s$coefficients))
  expect_equal(colnames(s$coefficients),
               c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  expect_equal(rownames(s$coefficients), c("mu", "log_sigma"))

  # Other components
  expect_equal(s$loglik, result$loglik)
  expect_true(is.finite(s$aic))
})


test_that("print methods work without error", {
  set.seed(42)
  x <- rnorm(50, mean = 3, sd = 1)

  result <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  # print.femtofit
  expect_output(print(result), "femtofit")
  expect_output(print(result), "Coefficients")

  # print.summary.femtofit
  expect_output(print(summary(result)), "femtofit")
  expect_output(print(summary(result)), "Log-likelihood")
  expect_output(print(summary(result)), "AIC")
})


test_that("observed_info() returns negative Hessian", {
  set.seed(42)
  x <- rnorm(50, mean = 3, sd = 1)

  result <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  info <- observed_info(result)
  expect_true(is.matrix(info))
  expect_equal(dim(info), c(2, 2))

  # Should be negative of Hessian
  expect_equal(info, -result$hessian)

  # Should be positive definite (at MLE of well-identified model)
  expect_true(all(eigen(info)$values > 0))
})


test_that("fit() works with Poisson distribution", {
  set.seed(42)
  x <- rpois(100, lambda = 3.5)

  result <- fit(
    function(lambda) loglik_poisson(lambda, x),
    params = c(lambda = 1)
  )

  expect_true(result$converged)
  expect_equal(unname(coef(result)["lambda"]), mean(x), tolerance = 0.1)
})


test_that("fit() works with exponential distribution", {
  set.seed(42)
  x <- rexp(100, rate = 2.5)

  result <- fit(
    function(rate) loglik_exponential(rate, x),
    params = c(rate = 1)
  )

  expect_true(result$converged)
  # MLE for exponential rate is 1/mean
  expect_equal(unname(coef(result)["rate"]), 1/mean(x), tolerance = 0.15)
})


test_that("fit() handles unnamed params by generating names", {
  set.seed(42)
  x <- rexp(50, rate = 2)

  result <- fit(
    function(p) loglik_exponential(p$p1, x),
    params = c(1)  # No names
  )

  expect_named(coef(result), "p1")
})


test_that("femtofit stores the original call", {
  set.seed(42)
  x <- rnorm(50, mean = 3, sd = 1)

  result <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  expect_false(is.null(result$call))
})


test_that("fit() produces reasonable standard errors", {
  set.seed(42)
  n <- 200
  true_mu <- 5
  true_sigma <- 2
  x <- rnorm(n, mean = true_mu, sd = true_sigma)

  result <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  # Theoretical SE for mu: sigma / sqrt(n)
  theoretical_se_mu <- true_sigma / sqrt(n)

  # Estimated SE should be in the right ballpark
  expect_equal(unname(se(result)["mu"]), theoretical_se_mu, tolerance = 0.05)

  # True value should be within 2 SEs of estimate (most of the time)
  expect_true(abs(coef(result)["mu"] - true_mu) < 2.5 * se(result)["mu"])
})
