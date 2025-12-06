test_that("lgamma works correctly", {
  x <- val(3)
  y <- lgamma(x)
  expect_equal(data(y), lgamma(3))

  backward(y)
  # d/dx lgamma(x) = digamma(x)
  expect_equal(grad(x), digamma(3))
})

test_that("digamma works correctly", {
  x <- val(2)
  y <- digamma(x)
  expect_equal(data(y), digamma(2))

  backward(y)
  # d/dx digamma(x) = trigamma(x)
  expect_equal(grad(x), trigamma(2))
})

test_that("trigamma works correctly", {
  x <- val(2.5)
  y <- trigamma(x)
  expect_equal(data(y), trigamma(2.5))

  backward(y)
  # d/dx trigamma(x) = psigamma(x, 2)
  expect_equal(grad(x), psigamma(2.5, 2))
})

test_that("log1p works correctly", {
  x <- val(0.5)
  y <- log1p(x)
  expect_equal(data(y), log1p(0.5))

  backward(y)
  # d/dx log(1+x) = 1/(1+x)
  expect_equal(grad(x), 1 / 1.5)
})

test_that("sin works correctly", {
  x <- val(pi / 4)
  y <- sin(x)
  expect_equal(data(y), sin(pi / 4))

  backward(y)
  # d/dx sin(x) = cos(x)
  expect_equal(grad(x), cos(pi / 4))
})

test_that("cos works correctly", {
  x <- val(pi / 3)
  y <- cos(x)
  expect_equal(data(y), cos(pi / 3))

  backward(y)
  # d/dx cos(x) = -sin(x)
  expect_equal(grad(x), -sin(pi / 3))
})

test_that("logit works correctly", {
  x <- val(0.3)
  y <- logit(x)
  expect_equal(data(y), log(0.3 / 0.7))

  backward(y)
  # d/dx logit(x) = 1/(x(1-x))
  expect_equal(grad(x), 1 / (0.3 * 0.7))
})

test_that("softplus works correctly", {
  x <- val(2)
  y <- softplus(x)
  expect_equal(data(y), log1p(exp(2)))

  backward(y)
  # d/dx softplus(x) = sigmoid(x)
  expected_grad <- 1 / (1 + exp(-2))
  expect_equal(grad(x), expected_grad)
})

test_that("chain rule works for composed functions", {
  # f(x) = sin(log(x)) at x = 2
  x <- val(2)
  y <- sin(log(x))

  backward(y)
  # df/dx = cos(log(x)) * (1/x)
  expected <- cos(log(2)) / 2
  expect_equal(grad(x), expected, tolerance = 1e-10)
})

test_that("lgamma gradient matches numerical gradient", {
  numerical_grad <- function(f, x, h = 1e-7) {
    (f(x + h) - f(x - h)) / (2 * h)
  }

  x_val <- 3.5
  x <- val(x_val)
  y <- lgamma(x)
  backward(y)

  numerical <- numerical_grad(lgamma, x_val)
  expect_equal(grad(x), numerical, tolerance = 1e-5)
})
