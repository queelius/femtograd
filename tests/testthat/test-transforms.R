# Tests for parameter transformation helpers

test_that("positive() transforms to positive values", {
  # With numeric
  expect_equal(positive(0), exp(0))
  expect_equal(positive(-2), exp(-2))
  expect_equal(positive(3), exp(3))

  # With value objects
  x <- val(-1)
  y <- positive(x)
  expect_equal(data(y), exp(-1))

  # Gradient flows correctly
  backward(y)
  # d/dx exp(x) = exp(x)
  expect_equal(grad(x), exp(-1))
})


test_that("probability() transforms to (0, 1)", {
  # With numeric
  expect_equal(probability(0), 0.5)
  expect_true(probability(-10) > 0 && probability(-10) < 0.5)
  expect_true(probability(10) > 0.5 && probability(10) < 1)

  # With value objects
  x <- val(2)
  p <- probability(x)
  expect_equal(data(p), sigmoid(2))

  # Gradient flows correctly
  backward(p)
  # d/dx sigmoid(x) = sigmoid(x) * (1 - sigmoid(x))
  s <- sigmoid(2)
  expect_equal(grad(x), s * (1 - s), tolerance = 1e-10)
})


test_that("bounded() transforms to (lower, upper)", {
  # With numeric
  expect_equal(bounded(0, 0, 10), 5)  # sigmoid(0) = 0.5
  expect_equal(bounded(0, -1, 1), 0)

  # Approaches bounds (use moderate values to avoid underflow)
  expect_true(bounded(-5, 0, 10) > 0 && bounded(-5, 0, 10) < 1)
  expect_true(bounded(5, 0, 10) > 9 && bounded(5, 0, 10) < 10)

  # With value objects
  x <- val(0)
  y <- bounded(x, 2, 8)
  expect_equal(data(y), 5)  # 2 + 6 * 0.5

  # Gradient flows correctly
  backward(y)
  # d/dx (a + (b-a)*sigmoid(x)) = (b-a) * sigmoid(x) * (1-sigmoid(x))
  s <- 0.5
  expect_equal(grad(x), 6 * s * (1 - s), tolerance = 1e-10)
})


test_that("lower_bounded() transforms to (lower, inf)", {
  # With numeric
  expect_equal(lower_bounded(0, 2), 3)  # 2 + exp(0) = 3
  expect_equal(lower_bounded(-1, 0), exp(-1))

  # Always above lower bound (use moderate value)
  expect_true(lower_bounded(-5, 5) > 5)

  # With value objects
  x <- val(1)
  y <- lower_bounded(x, 3)
  expect_equal(data(y), 3 + exp(1))

  backward(y)
  expect_equal(grad(x), exp(1), tolerance = 1e-10)
})


test_that("upper_bounded() transforms to (-inf, upper)", {
  # With numeric
  expect_equal(upper_bounded(0, 10), 9)  # 10 - exp(0) = 9
  expect_equal(upper_bounded(1, 5), 5 - exp(1))

  # Always below upper bound (use moderate value)
  expect_true(upper_bounded(-5, 5) < 5)

  # With value objects
  x <- val(0)
  y <- upper_bounded(x, 7)
  expect_equal(data(y), 6)  # 7 - 1

  backward(y)
  # d/dx (upper - exp(x)) = -exp(x)
  expect_equal(grad(x), -1, tolerance = 1e-10)
})


test_that("inverse transforms recover original values", {
  # inv_positive
  expect_equal(inv_positive(exp(3)), 3)
  expect_equal(inv_positive(1), 0)

  # inv_probability
  expect_equal(inv_probability(0.5), 0, tolerance = 1e-10)
  expect_equal(inv_probability(sigmoid(2)), 2, tolerance = 1e-10)

  # inv_bounded
  expect_equal(inv_bounded(5, 0, 10), 0, tolerance = 1e-10)  # 5 is midpoint
  expect_equal(inv_bounded(0, -1, 1), 0, tolerance = 1e-10)
})


test_that("transforms compose correctly in optimization", {
  set.seed(42)
  # Generate data from known parameters
  true_mu <- 5
  true_sigma <- 2
  x <- rnorm(100, mean = true_mu, sd = true_sigma)

  # Fit using positive() for sigma constraint
  result <- fit(
    function(mu, log_sigma) {
      sigma <- positive(log_sigma)
      loglik_normal(mu, sigma, x)
    },
    params = c(mu = 0, log_sigma = 0)
  )

  expect_true(result$converged)

  # Recover sigma on natural scale
  fitted_sigma <- positive(coef(result)["log_sigma"])
  expect_equal(unname(fitted_sigma), true_sigma, tolerance = 0.3)
})


test_that("bounded() validates inputs", {
  expect_error(bounded(0, 10, 5), "lower < upper")  # Invalid bounds
})


test_that("transforms work with vectors", {
  # positive with vector
  x <- c(-1, 0, 1)
  expect_equal(positive(x), exp(x))

  # probability with vector
  expect_equal(probability(x), sigmoid(x))

  # bounded with vector
  expect_equal(bounded(x, 0, 10), 10 * sigmoid(x))
})


test_that("inverse transforms are true inverses", {
  # Test round-trip for various values
  test_vals <- c(-3, -1, 0, 1, 3)

  for (v in test_vals) {
    # positive <-> inv_positive
    expect_equal(inv_positive(positive(v)), v, tolerance = 1e-10)

    # probability <-> inv_probability
    expect_equal(inv_probability(probability(v)), v, tolerance = 1e-10)

    # bounded <-> inv_bounded
    expect_equal(inv_bounded(bounded(v, 2, 8), 2, 8), v, tolerance = 1e-10)
  }
})
