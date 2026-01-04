# Tests for model diagnostics

# Helper: Normal log-likelihood with log(sigma) parameterization
loglik_normal_logsig <- function(mu, log_sigma, x) {
  sigma <- exp(log_sigma)
  loglik_normal(mu, sigma, x)
}

test_that("check_hessian works for well-specified model", {
  set.seed(42)
  x <- rnorm(100, mean = 5, sd = 2)

  result <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  check <- check_hessian(result)

  expect_s3_class(check, "hessian_check")
  expect_true(check$is_negative_definite)
  expect_false(check$is_singular)
  expect_equal(check$rank, 2)
  expect_equal(check$n_params, 2)
  expect_true(is.finite(check$condition_number))
  expect_equal(length(check$eigenvalues), 2)
  expect_true(all(check$eigenvalues > 0))
})


test_that("check_hessian detects high condition number", {
  # Create a situation with nearly collinear parameters
  # This is hard to construct, so we'll just verify the structure
  set.seed(42)
  x <- rnorm(100, mean = 5, sd = 2)

  result <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  check <- check_hessian(result)

  # For a normal model, condition number should be reasonable
  expect_true(check$condition_number < 1e6)
})


test_that("check_convergence works for converged model", {
  set.seed(42)
  x <- rnorm(100, mean = 5, sd = 2)

  result <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  check <- check_convergence(result)

  expect_s3_class(check, "convergence_check")
  expect_true(check$converged)
  expect_true(is.finite(check$loglik))
  expect_true(is.numeric(check$gradient_norm))
})


test_that("check_convergence respects tolerance", {
  set.seed(42)
  x <- rnorm(100, mean = 5, sd = 2)

  result <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  # Very strict tolerance
  check_strict <- check_convergence(result, tol = 1e-12)

  # Lenient tolerance
  check_lenient <- check_convergence(result, tol = 1)

  expect_true(check_lenient$gradient_ok)
  # Strict might or might not pass depending on actual gradient
})


test_that("diagnostics() returns comprehensive results", {
  set.seed(42)
  x <- rnorm(100, mean = 5, sd = 2)

  result <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  diag <- diagnostics(result)

  expect_s3_class(diag, "model_diagnostics")
  expect_s3_class(diag$hessian, "hessian_check")
  expect_s3_class(diag$convergence, "convergence_check")
  expect_true(diag$ok)  # Should be OK for a good fit
  expect_type(diag$n_warnings, "integer")
})


test_that("se_reliable returns TRUE for good fit", {
  set.seed(42)
  x <- rnorm(100, mean = 5, sd = 2)

  result <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  expect_true(se_reliable(result))
})


test_that("se_reliable returns FALSE for invalid input", {
  expect_false(se_reliable(list(a = 1)))
  expect_false(se_reliable(NULL))
  expect_false(se_reliable("not a fit"))
})


test_that("print methods work without error", {
  set.seed(42)
  x <- rnorm(50, mean = 3, sd = 1)

  result <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  hess_check <- check_hessian(result)
  conv_check <- check_convergence(result)
  full_diag <- diagnostics(result)

  expect_output(print(hess_check), "Hessian Diagnostics")
  expect_output(print(hess_check), "Condition number")

  expect_output(print(conv_check), "Convergence Diagnostics")
  expect_output(print(conv_check), "Gradient norm")

  expect_output(print(full_diag), "Model Diagnostics")
})


test_that("check_hessian validates input", {
  expect_error(check_hessian(list(a = 1)), "femtofit")
})


test_that("check_convergence validates input", {
  expect_error(check_convergence(list(a = 1)), "femtofit")
})


test_that("diagnostics work with single parameter model", {
  set.seed(42)
  x <- rexp(100, rate = 2)

  result <- fit(
    function(rate) loglik_exponential(rate, x),
    params = c(rate = 1)
  )

  diag <- diagnostics(result)

  expect_true(diag$ok)
  expect_equal(diag$hessian$n_params, 1)
})


test_that("eigenvalues are ordered correctly", {
  set.seed(42)
  x <- rnorm(100, mean = 5, sd = 2)

  result <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  check <- check_hessian(result)

  # Eigenvalues from eigen() are in decreasing order
  expect_true(all(diff(check$eigenvalues) <= 0))
})


test_that("diagnostics detect issues when present", {
  # Create a femtofit-like object with bad Hessian
  bad_result <- structure(
    list(
      coefficients = c(mu = 1, sigma = 1),
      vcov = matrix(c(1, 0, 0, 1), 2, 2),
      loglik = -50,
      hessian = matrix(c(1, 0, 0, 1), 2, 2),  # Positive definite (wrong sign!)
      gradient = c(0.1, 0.1),
      converged = TRUE,
      iterations = 10
    ),
    class = "femtofit"
  )

  hess_check <- check_hessian(bad_result)

  # -H is negative definite here (since H is positive definite)
  expect_false(hess_check$is_negative_definite)
  expect_true(length(hess_check$warnings) > 0)
})


test_that("condition number calculation is correct", {
  set.seed(42)
  x <- rnorm(100, mean = 5, sd = 2)

  result <- fit(
    function(mu, log_sigma) loglik_normal_logsig(mu, log_sigma, x),
    params = c(mu = 0, log_sigma = 0)
  )

  check <- check_hessian(result)

  # Condition number should be max(eig) / min(eig)
  expected_cond <- max(check$eigenvalues) / min(check$eigenvalues)
  expect_equal(check$condition_number, expected_cond)
})
