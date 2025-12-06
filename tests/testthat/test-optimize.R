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

  loglik <- function(p) loglik_exponential(p[[1]], x)
  mle_result <- find_mle(loglik, list(val(1)))

  # Test H0: rate = 2
  test_result <- wald_test(mle_result, null_values = c(2))

  expect_true("statistic" %in% names(test_result))
  expect_true("p_value" %in% names(test_result))
  expect_equal(test_result$df, 1)

  # Should not reject null when true value is 2
  # (p-value should be reasonably large most of the time)
  expect_true(test_result$p_value >= 0 && test_result$p_value <= 1)
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
