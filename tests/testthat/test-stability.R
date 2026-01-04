test_that("logsumexp handles normal values correctly", {
  # Simple case
  x <- c(1, 2, 3)
  result <- logsumexp(x)
  expected <- log(sum(exp(x)))
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("logsumexp prevents overflow with large values", {
  # exp(1000) would overflow
  x <- c(1000, 1001, 1002)
  result <- logsumexp(x)

  # Should be approximately 1002 + log(1 + exp(-1) + exp(-2))
  expected <- 1002 + log(1 + exp(-1) + exp(-2))
  expect_equal(result, expected, tolerance = 1e-10)
  expect_true(is.finite(result))
})

test_that("logsumexp handles very negative values", {
  x <- c(-1000, -1001, -1002)
  result <- logsumexp(x)

  # Should be approximately -1000 + log(1 + exp(-1) + exp(-2))
  expected <- -1000 + log(1 + exp(-1) + exp(-2))
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("logsumexp works with value objects", {
  x <- val(c(1, 2, 3))
  result <- logsumexp(x)

  expected <- log(sum(exp(c(1, 2, 3))))
  expect_equal(data(result), expected, tolerance = 1e-10)

  # Check gradients
  backward(result)
  # d/dx_i logsumexp(x) = exp(x_i) / sum(exp(x)) = softmax(x)_i
  expected_grad <- exp(c(1, 2, 3)) / sum(exp(c(1, 2, 3)))
  expect_equal(as.vector(grad(x)), expected_grad, tolerance = 1e-10)
})

test_that("logsumexp handles edge cases", {
  # Single value
  expect_equal(logsumexp(5), 5)

  # All -Inf
  expect_equal(logsumexp(c(-Inf, -Inf)), -Inf)

  # Empty (numeric version)
  expect_equal(logsumexp(numeric(0)), -Inf)
})

test_that("softmax produces valid probabilities", {
  x <- c(1, 2, 3)
  result <- softmax(x)

  expect_equal(sum(result), 1, tolerance = 1e-10)
  expect_true(all(result > 0))
  expect_true(all(result < 1))
})

test_that("softmax is stable for large values", {
  x <- c(1000, 1001, 1002)
  result <- softmax(x)

  expect_equal(sum(result), 1, tolerance = 1e-10)
  expect_true(all(is.finite(result)))
})

test_that("softmax works with value objects", {
  x <- val(c(1, 2, 3))
  result <- softmax(x)

  expect_equal(sum(data(result)), 1, tolerance = 1e-10)
})

test_that("log_safe handles zeros", {
  x <- c(0, 1, 2)
  result <- log_safe(x)

  expect_true(is.finite(result[1]))  # log(0) would be -Inf
  expect_equal(result[2], 0)
  expect_equal(result[3], log(2))
})

test_that("log_safe works with value objects", {
  x <- val(c(0.001, 1, 2))
  result <- log_safe(x)

  expect_true(all(is.finite(data(result))))

  backward(result)
  expect_true(all(is.finite(grad(x))))
})

test_that("div_safe handles division by zero", {
  result <- div_safe(1, 0)
  expect_true(is.finite(result))

  result <- div_safe(c(1, 2, 3), c(0, 1, 0))
  expect_true(all(is.finite(result)))
})

test_that("div_safe works with value objects", {
  x <- val(1)
  y <- val(0.0001)
  result <- div_safe(x, y)

  expect_true(is.finite(data(result)))

  backward(result)
  expect_true(is.finite(grad(x)))
  expect_true(is.finite(grad(y)))
})

test_that("sigmoid_stable matches sigmoid for normal values", {
  x <- c(-2, -1, 0, 1, 2)
  result <- sigmoid_stable(x)
  expected <- 1 / (1 + exp(-x))

  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("sigmoid_stable handles extreme values", {
  # Large positive
  result_pos <- sigmoid_stable(1000)
  expect_equal(result_pos, 1, tolerance = 1e-10)

  # Large negative
  result_neg <- sigmoid_stable(-1000)
  expect_equal(result_neg, 0, tolerance = 1e-10)

  # Both should be finite
  expect_true(is.finite(result_pos))
  expect_true(is.finite(result_neg))
})

test_that("sigmoid_stable works with value objects", {
  x <- val(0)
  result <- sigmoid_stable(x)

  expect_equal(data(result), 0.5)

  backward(result)
  expect_equal(grad(x), 0.25)  # sigmoid'(0) = 0.5 * 0.5
})

test_that("exp_safe prevents overflow", {
  # exp(710) would overflow
  result <- exp_safe(1000)
  expect_true(is.finite(result))
  expect_equal(result, exp(709))  # Clamped to max_val
})

test_that("exp_safe works with value objects", {
  x <- val(2)
  result <- exp_safe(x)

  expect_equal(data(result), exp(2))

  backward(result)
  expect_equal(grad(x), exp(2))
})

test_that("log_sigmoid is stable for extreme values", {
  # Large positive: log(sigmoid(1000)) should be close to 0
  result_pos <- log_sigmoid(1000)
  expect_equal(result_pos, 0, tolerance = 1e-10)

  # Large negative: log(sigmoid(-1000)) should be close to -1000
  result_neg <- log_sigmoid(-1000)
  expect_equal(result_neg, -1000, tolerance = 1)

  expect_true(is.finite(result_pos))
  expect_true(is.finite(result_neg))
})

test_that("log_sigmoid works with value objects", {
  x <- val(0)
  result <- log_sigmoid(x)

  # log(sigmoid(0)) = log(0.5)
  expect_equal(data(result), log(0.5), tolerance = 1e-10)

  backward(result)
  # d/dx log(sigmoid(x)) = 1 - sigmoid(x) = 0.5 at x=0
  expect_equal(grad(x), 0.5, tolerance = 1e-10)
})

test_that("log_sigmoid matches direct computation for moderate values", {
  x <- c(-5, -2, 0, 2, 5)

  result <- log_sigmoid(x)
  expected <- log(1 / (1 + exp(-x)))

  expect_equal(result, expected, tolerance = 1e-10)
})
