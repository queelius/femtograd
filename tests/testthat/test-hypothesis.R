# Tests for hypothesis testing functions

# Helper: Normal log-likelihood with log(sigma) parameterization
loglik_normal_logsig <- function(mu, log_sigma, x) {
  sigma <- exp(log_sigma)
  loglik_normal(mu, sigma, x)
}

test_that("lrt() works with femtofit objects", {
  set.seed(42)
  x <- rnorm(100, mean = 5, sd = 2)

  # Full model: estimate mu and sigma
  full <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  # Null model: mu fixed at 0 (wrong model)
  null <- fit(
    function(log_sigma) loglik_normal_logsig(0, log_sigma, x),
    params = c(log_sigma = 0)
  )

  test <- lrt(null, full)

  # Check structure
  expect_s3_class(test, "likelihood_ratio_test")
  expect_s3_class(test, "hypothesis_test")
  expect_true(is.numeric(test$stat))
  expect_true(is.numeric(test$p.value))
  expect_equal(test$dof, 1)  # 1 extra parameter in full model

  # Statistic should be positive
  expect_true(test$stat >= 0)

  # P-value in [0, 1]
  expect_true(test$p.value >= 0 && test$p.value <= 1)

  # Since true mu = 5 != 0, the test should be significant
  expect_true(is_significant_at(test, 0.05))
})


test_that("lrt() works with raw log-likelihoods", {
  test <- lrt(null_model = -150, alt_model = -140, df = 2)

  expect_s3_class(test, "likelihood_ratio_test")
  expect_equal(test$stat, 20)  # -2 * (-150 - (-140)) = 20
  expect_equal(test$dof, 2)
  expect_equal(test$null_loglik, -150)
  expect_equal(test$alt_loglik, -140)
})


test_that("lrt() validates inputs", {
  # Missing df with raw values
  expect_error(lrt(-150, -140), "df must be specified")

  # Wrong model order (null has more params)
  set.seed(42)
  x <- rnorm(50)

  simple <- fit(
    function(log_sigma) loglik_normal_logsig(0, log_sigma, x),
    params = c(log_sigma = 0)
  )

  complex <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  expect_error(lrt(complex, simple), "more parameters")
})


test_that("wald_test.femtofit() tests single parameter", {
  set.seed(42)
  x <- rnorm(100, mean = 5, sd = 2)

  result <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  # Test mu = 0 (should reject since true mu = 5)
  test <- wald_test(result, "mu")

  expect_s3_class(test, "wald_test")
  expect_s3_class(test, "hypothesis_test")
  expect_equal(test$dof, 1)
  expect_true(is_significant_at(test, 0.05))

  # Test mu = 5 (should not reject)
  test_correct <- wald_test(result, "mu", null_value = 5)
  expect_false(is_significant_at(test_correct, 0.05))
})


test_that("wald_test.femtofit() tests all parameters", {
  set.seed(42)
  x <- rnorm(50, mean = 3, sd = 1)

  result <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  # Test all parameters
  tests <- wald_test(result)

  expect_true(is.list(tests))
  expect_equal(length(tests), 2)
  expect_named(tests, c("mu", "log_sigma"))

  expect_s3_class(tests$mu, "wald_test")
  expect_s3_class(tests$log_sigma, "wald_test")
})


test_that("wald_test.default() works with raw values", {
  test <- wald_test(2.5, se = 0.8, null_value = 0)

  expect_s3_class(test, "wald_test")
  expect_equal(test$estimate, 2.5)
  expect_equal(test$se, 0.8)
  expect_equal(test$null_value, 0)
  expect_equal(test$z, 2.5 / 0.8)
  expect_equal(test$stat, (2.5 / 0.8)^2)
})


test_that("pval() extracts p-value correctly", {
  test <- wald_test(2.5, se = 0.8)
  expect_equal(pval(test), test$p.value)

  lrt_test <- lrt(-150, -140, df = 1)
  expect_equal(pval(lrt_test), lrt_test$p.value)
})


test_that("test_stat() extracts statistic correctly", {
  test <- wald_test(2.5, se = 0.8)
  expect_equal(test_stat(test), test$stat)
})


test_that("dof() extracts degrees of freedom correctly", {
  test <- wald_test(2.5, se = 0.8)
  expect_equal(dof(test), 1)

  lrt_test <- lrt(-150, -140, df = 3)
  expect_equal(dof(lrt_test), 3)
})


test_that("is_significant_at() works correctly", {
  # Very significant test
  significant <- wald_test(5.0, se = 0.5)  # z = 10
  expect_true(is_significant_at(significant, 0.05))
  expect_true(is_significant_at(significant, 0.01))
  expect_true(is_significant_at(significant, 0.001))

  # Not significant test
  not_sig <- wald_test(0.1, se = 0.5)  # z = 0.2
  expect_false(is_significant_at(not_sig, 0.05))
})


test_that("print methods work without error", {
  test_wald <- wald_test(2.5, se = 0.8)
  test_lrt <- lrt(-150, -140, df = 2)

  expect_output(print(test_wald), "Wald Test")
  expect_output(print(test_wald), "z-score")

  expect_output(print(test_lrt), "Likelihood Ratio Test")
  expect_output(print(test_lrt), "LRT statistic")
})


test_that("hypothesize compatibility: structure matches", {
  # Tests should have the fields expected by hypothesize package
  test <- wald_test(2.5, se = 0.8)

  expect_true("stat" %in% names(test))
  expect_true("p.value" %in% names(test))
  expect_true("dof" %in% names(test))
  expect_true("hypothesis_test" %in% class(test))
})


test_that("lrt() gives sensible results for nested models", {
  set.seed(123)

  # Generate data from exponential(rate=2)
  x <- rexp(200, rate = 2)

  # Full model: estimate rate
  full <- fit(
    function(rate) loglik_exponential(rate, x),
    params = c(rate = 1)
  )

  # Null model: rate = 1 (wrong)
  # We compute null loglik directly
  null_ll <- sum(-x)  # loglik_exponential with rate=1

  test <- lrt(null_model = null_ll, alt_model = full$loglik, df = 1)

  # True rate is 2, so rate=1 should be rejected
  expect_true(is_significant_at(test, 0.05))
})


test_that("wald_test results match summary z-values", {
  set.seed(42)
  x <- rnorm(100, mean = 5, sd = 2)

  result <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  summ <- summary(result)
  tests <- wald_test(result)

  # z-values should match
  expect_equal(tests$mu$z, summ$coefficients["mu", "z value"],
               tolerance = 1e-10)
  expect_equal(tests$log_sigma$z, summ$coefficients["log_sigma", "z value"],
               tolerance = 1e-10)
})
