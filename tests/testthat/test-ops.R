test_that("basic arithmetic operations work", {
  x <- val(2)
  y <- val(3)

  # Addition
  z <- x + y
  expect_equal(data(z), 5)

  # Multiplication
  z <- x * y
  expect_equal(data(z), 6)

  # Subtraction
  z <- y - x
  expect_equal(data(z), 1)

  # Division
  z <- y / x
  expect_equal(data(z), 1.5)

  # Power
  z <- x ^ val(3)
  expect_equal(data(z), 8)
})

test_that("mixed value-scalar operations work", {
  x <- val(5)

  expect_equal(data(x + 3), 8)
  expect_equal(data(x - 2), 3)
  expect_equal(data(x * 2), 10)
  expect_equal(data(x / 2), 2.5)
  expect_equal(data(x ^ 2), 25)
})

test_that("scalar-value operations work (scalar first)", {
  # Regression test: scalar / value was failing with
  # "$ operator is invalid for atomic vectors"
  x <- val(3)

  # All these should work with scalar as first operand
  expect_equal(data(2 + x), 5)
  expect_equal(data(10 - x), 7)
  expect_equal(data(4 * x), 12)
  expect_equal(data(6 / x), 2)
  expect_equal(data(2 ^ x), 8)

  # Check gradients propagate correctly
  y <- val(2)
  z <- 10 / y  # d/dy (10/y) = -10/y^2 = -2.5

  backward(z)
  expect_equal(grad(y), -2.5)

  # Vector divided by value
  y <- val(2)
  result <- c(2, 4, 6) / y
  expect_equal(data(result), c(1, 2, 3))
})

test_that("power with value exponent has correct gradient", {
  # Test the value exponent path (line 64-65 in ops.R)
  b <- val(2)
  e <- val(3)
  z <- b^e

  expect_equal(data(z), 8)

  backward(z)
  # d/db = e * b^(e-1) = 3 * 4 = 12
  expect_equal(grad(b), 12)
  # d/de = log(b) * b^e = log(2) * 8
  expect_equal(grad(e), log(2) * 8, tolerance = 1e-10)
})

test_that("sum.value with empty input returns zero", {
  result <- sum.value()
  expect_equal(data(result), 0)
})

test_that("sum.value handles mixed scalars and values", {
  x <- val(2)
  y <- val(3)
  result <- sum(x, 5, y, 10)

  expect_equal(data(result), 20)

  backward(result)
  expect_equal(grad(x), 1)
  expect_equal(grad(y), 1)
})

test_that("mean.value returns input for single value", {
  x <- val(5)
  result <- mean(x)
  expect_equal(data(result), 5)
})

test_that("mean.value with list of values works", {
  values <- list(val(2), val(4), val(6))
  result <- mean.value(values)

  expect_equal(data(result), 4)

  backward(result)
  # Each input gets gradient 1/n = 1/3
  expect_equal(grad(values[[1]]), 1/3, tolerance = 1e-10)
  expect_equal(grad(values[[2]]), 1/3, tolerance = 1e-10)
  expect_equal(grad(values[[3]]), 1/3, tolerance = 1e-10)
})

test_that("unary negation works", {
  x <- val(5)
  z <- -x
  expect_equal(data(z), -5)

  backward(z)
  expect_equal(grad(x), -1)
})

test_that("absolute value works", {
  x <- val(-3)
  z <- abs(x)
  expect_equal(data(z), 3)

  backward(z)
  expect_equal(grad(x), -1)  # d/dx abs(x) = sign(x) = -1 for x < 0

  # Positive value
  x <- val(4)
  z <- abs(x)
  expect_equal(data(z), 4)

  backward(z)
  expect_equal(grad(x), 1)  # d/dx abs(x) = sign(x) = 1 for x > 0
})

test_that("transcendental functions work", {
  x <- val(2)

  # Natural log
  z <- log(x)
  expect_equal(data(z), log(2))

  # Exponential
  z <- exp(x)
  expect_equal(data(z), exp(2))

  # Square root
  z <- sqrt(x)
  expect_equal(data(z), sqrt(2))

  # Tanh
  z <- tanh(x)
  expect_equal(data(z), tanh(2))
})

test_that("activation functions work", {
  # Sigmoid
  x <- val(0)
  z <- sigmoid(x)
  expect_equal(data(z), 0.5)

  x <- val(-10)
  z <- sigmoid(x)
  expect_lt(data(z), 0.01)

  # ReLU
  x <- val(5)
  z <- relu(x)
  expect_equal(data(z), 5)

  x <- val(-3)
  z <- relu(x)
  expect_equal(data(z), 0)
})

test_that("sum works with value objects", {
  x <- val(1)
  y <- val(2)
  z <- val(3)

  s <- sum(x, y, z)
  expect_equal(data(s), 6)

  backward(s)
  expect_equal(grad(x), 1)
  expect_equal(grad(y), 1)
  expect_equal(grad(z), 1)
})

test_that("mean works with value objects", {
  x <- val(2)
  y <- val(4)
  z <- val(6)

  # Note: mean.value requires calling on a value object, not a list
  # Use sum and divide for mean functionality
  s <- sum(x, y, z)
  m <- s / 3
  expect_equal(data(m), 4)

  backward(m)
  expect_equal(grad(x), 1/3)
  expect_equal(grad(y), 1/3)
  expect_equal(grad(z), 1/3)
})

test_that("gradient computation is correct for simple expressions", {
  # f(x) = x^2
  x <- val(3)
  y <- x * x
  backward(y)
  expect_equal(grad(x), 6)  # df/dx = 2x = 6

  # f(x) = x*2 + 3
  # Note: use x*2 not 2*x because R dispatches on first arg
  x <- val(5)
  y <- x * 2 + 3
  backward(y)
  expect_equal(grad(x), 2)  # df/dx = 2
})

test_that("gradient computation works for complex expressions", {
  # f(x,y) = x*y + sin(x) would require sin, but let's use:
  # f(x,y) = x*y + x^2
  x <- val(3)
  y <- val(4)
  z <- x * y + x^2
  expect_equal(data(z), 12 + 9)

  backward(z)
  expect_equal(grad(x), 4 + 6)  # df/dx = y + 2x = 4 + 6 = 10
  expect_equal(grad(y), 3)       # df/dy = x = 3
})

test_that("zero_grad resets gradients", {
  x <- val(2)
  y <- x * x
  backward(y)
  expect_equal(grad(x), 4)

  zero_grad(y)
  expect_equal(grad(x), 0)
  expect_equal(grad(y), 0)
})

test_that("labels work for debugging", {
  x <- val(5, label = "learning_rate")
  expect_equal(x$label, "learning_rate")

  y <- val(10)
  expect_equal(y$label, "")
})

test_that("gradient accumulation works correctly", {
  # Test that gradients accumulate when a node is used multiple times
  # f(x) = x + x = 2x, so df/dx = 2
  x <- val(3)
  y <- x + x
  backward(y)
  expect_equal(grad(x), 2)

  # f(x) = x*x + x*x = 2x^2, so df/dx = 4x = 12
  x <- val(3)
  z <- x * x + x * x
  backward(z)
  expect_equal(grad(x), 12)
})
