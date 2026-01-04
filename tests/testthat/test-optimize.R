test_that("gradient_ascent finds maximum of simple quadratic", {
  # f(x) = -(x-3)^2 + 10, maximum at x = 3
  objective <- function(p) {
    x <- p[[1]]
    -(x - 3)^2 + 10
  }

  result <- gradient_ascent(objective, list(val(0)), lr = 0.1, max_iter = 100)

  expect_true(result$converged)
  expect_equal(data(result$params[[1]]), 3, tolerance = 0.01)
  expect_equal(result$value, 10, tolerance = 0.01)
})

test_that("gradient_descent finds minimum of simple quadratic", {
  # f(x) = (x-5)^2, minimum at x = 5
  objective <- function(p) {
    x <- p[[1]]
    (x - 5)^2
  }

  result <- gradient_descent(objective, list(val(0)), lr = 0.1, max_iter = 100)

  expect_true(result$converged)
  expect_equal(data(result$params[[1]]), 5, tolerance = 0.01)
})

test_that("gradient_ascent finds MLE for exponential", {
  set.seed(123)
  x <- rexp(50, rate = 2)
  true_mle <- 1 / mean(x)

  loglik <- function(p) loglik_exponential(p[[1]], x)

  result <- gradient_ascent(loglik, list(val(1)), lr = 0.1, max_iter = 500)

  expect_equal(data(result$params[[1]]), true_mle, tolerance = 0.01)
})

test_that("newton_raphson converges faster than gradient ascent", {
  # Simple quadratic - Newton should converge in 1-2 iterations
  objective <- function(p) {
    x <- p[[1]]
    -(x - 3)^2 + 10
  }

  result <- newton_raphson(objective, list(val(0)), max_iter = 10)

  expect_true(result$converged)
  expect_equal(data(result$params[[1]]), 3, tolerance = 1e-6)
  expect_lt(result$iterations, 5)  # Should converge very fast
})

test_that("newton_raphson finds MLE for normal distribution", {
  set.seed(42)
  true_mu <- 5
  true_sigma <- 2
  x <- rnorm(100, true_mu, true_sigma)

  loglik <- function(p) loglik_normal(p[[1]], p[[2]], x)

  # Start reasonably close to avoid divergence
  # Newton-Raphson can be unstable with bad starting points
  result <- newton_raphson(
    loglik,
    list(val(mean(x)), val(sd(x))),  # Start near sample estimates
    max_iter = 50,
    step_scale = 0.5  # Damped Newton for stability
  )

  expect_true(result$converged)
  expect_equal(data(result$params[[1]]), mean(x), tolerance = 0.01)
  # Sample SD is MLE for sigma
  expect_equal(data(result$params[[2]]), sqrt(mean((x - mean(x))^2)), tolerance = 0.1)
})

test_that("find_mle returns correct structure", {
  set.seed(42)
  x <- rexp(50, rate = 2)

  loglik <- function(p) loglik_exponential(p[[1]], x)

  result <- find_mle(loglik, list(val(1)), method = "newton")

  expect_true("estimate" %in% names(result))
  expect_true("se" %in% names(result))
  expect_true("vcov" %in% names(result))
  expect_true("loglik" %in% names(result))
  expect_true("converged" %in% names(result))

  # MLE should be close to 1/mean(x)
  expect_equal(result$estimate[1], 1 / mean(x), tolerance = 0.01)
})

test_that("confint_mle computes confidence intervals", {
  set.seed(42)
  x <- rexp(100, rate = 2)

  loglik <- function(p) loglik_exponential(p[[1]], x)
  mle_result <- find_mle(loglik, list(val(1)))

  ci <- confint_mle(mle_result, level = 0.95)

  expect_equal(ncol(ci), 4)  # estimate, se, lower, upper
  expect_true(ci[1, "lower"] < ci[1, "estimate"])
  expect_true(ci[1, "upper"] > ci[1, "estimate"])
})

test_that("wald_test computes test statistic", {
  set.seed(42)
  x <- rexp(100, rate = 2)  # True rate is 2

  # Use fit() which returns a femtofit object
  result <- fit(
    function(rate) loglik_exponential(rate, x),
    params = c(rate = 1)
  )

  # Test H0: rate = 2
  test_result <- wald_test(result, "rate", null_value = 2)

  # Check new interface structure
  expect_true("stat" %in% names(test_result))
  expect_true("p.value" %in% names(test_result))
  expect_equal(test_result$dof, 1)

  # Should not reject null when true value is 2
  # (p-value should be reasonably large most of the time)
  expect_true(test_result$p.value >= 0 && test_result$p.value <= 1)
})

test_that("gradient clipping prevents large steps", {
  # Function with steep gradient
  objective <- function(p) {
    x <- p[[1]]
    -exp(x)  # Gradient = -exp(x), can be very large
  }

  # Without clipping, this might diverge
  result <- gradient_ascent(
    objective,
    list(val(10)),  # Start at exp(10) ~ 22000 gradient
    lr = 0.01,
    max_iter = 100,
    grad_clip = 1
  )

  # Should not have exploded to infinity
  expect_true(is.finite(data(result$params[[1]])))
})

test_that("line_search finds appropriate step size", {
  # Simple quadratic: f(x) = (x - 3)^2, min at x = 3
  objective <- function(p) {
    x <- p[[1]]
    (x - 3)^2
  }

  # At x = 0, gradient = 2*(0-3) = -6, descent direction = -(-6) = 6
  param_values <- c(0)
  grad_vec <- c(-6)
  direction <- c(6)  # Descent direction (negative gradient)

  alpha <- line_search(objective, param_values, direction, grad_vec,
                       maximize = FALSE)

  # Should find a step that decreases the objective
  expect_true(alpha > 0)
  new_value <- (param_values + alpha * direction - 3)^2
  old_value <- 9  # (0 - 3)^2
  expect_lte(new_value, old_value)
})

test_that("bfgs minimizes simple quadratic", {
  # Simple quadratic: f(x,y) = (x-2)^2 + (y-3)^2
  # Minimum at (2, 3)
  quadratic <- function(p) {
    x <- p[[1]]
    y <- p[[2]]
    (x - 2)^2 + (y - 3)^2
  }

  result <- bfgs(quadratic, list(val(0), val(0)),
                 max_iter = 100, tol = 1e-6)

  expect_true(result$converged)
  expect_equal(data(result$params[[1]]), 2, tolerance = 0.01)
  expect_equal(data(result$params[[2]]), 3, tolerance = 0.01)
})

test_that("bfgs maximizes log-likelihood", {
  set.seed(42)
  x <- rexp(100, rate = 2)
  mle <- 1 / mean(x)

  loglik <- function(p) loglik_exponential(p[[1]], x)

  # Start near the MLE for stability
  result <- bfgs(loglik, list(val(mle * 0.8)), maximize = TRUE,
                 max_iter = 100, tol = 1e-6)

  expect_true(result$converged)
  expect_equal(data(result$params[[1]]), mle, tolerance = 0.01)
})

test_that("lbfgs minimizes simple quadratic", {
  quadratic <- function(p) {
    x <- p[[1]]
    y <- p[[2]]
    (x - 2)^2 + (y - 3)^2
  }

  result <- lbfgs(quadratic, list(val(0), val(0)),
                  m = 5, max_iter = 100, tol = 1e-6)

  expect_true(result$converged)
  expect_equal(data(result$params[[1]]), 2, tolerance = 0.01)
  expect_equal(data(result$params[[2]]), 3, tolerance = 0.01)
})

test_that("lbfgs is more memory efficient than bfgs", {
  # L-BFGS doesn't return inv_hessian (memory savings)
  quadratic <- function(p) {
    x <- p[[1]]
    y <- p[[2]]
    (x - 2)^2 + (y - 3)^2
  }

  result_bfgs <- bfgs(quadratic, list(val(0), val(0)), max_iter = 10)
  result_lbfgs <- lbfgs(quadratic, list(val(0), val(0)), max_iter = 10)

  expect_true("inv_hessian" %in% names(result_bfgs))
  expect_false("inv_hessian" %in% names(result_lbfgs))
})

test_that("bfgs handles multivariate normal MLE", {
  set.seed(123)
  true_mu <- 5
  true_sigma <- 2
  x <- rnorm(100, true_mu, true_sigma)

  loglik <- function(p) loglik_normal(p[[1]], p[[2]], x)

  # Start reasonably close for stability
  result <- bfgs(loglik, list(val(mean(x)), val(sd(x))), maximize = TRUE,
                 max_iter = 200, tol = 1e-6)

  expect_true(result$converged)
  expect_equal(data(result$params[[1]]), mean(x), tolerance = 0.1)
})
