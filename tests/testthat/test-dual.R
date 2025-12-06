test_that("dual number creation works", {
  d <- dual_num(3, 1)
  expect_true(is_dual(d))
  expect_equal(primal(d), 3)
  expect_equal(tangent(d), 1)
})

test_that("primal/tangent extract correctly from non-dual", {
  expect_equal(primal(5), 5)
  expect_equal(tangent(5), 0)
})

test_that("dual addition works", {
  d1 <- dual_num(2, 1)
  d2 <- dual_num(3, 2)
  result <- d1 + d2
  expect_equal(primal(result), 5)
  expect_equal(tangent(result), 3)
})

test_that("dual subtraction works", {
  d1 <- dual_num(5, 3)
  d2 <- dual_num(2, 1)
  result <- d1 - d2
  expect_equal(primal(result), 3)
  expect_equal(tangent(result), 2)
})

test_that("dual unary negation works", {
  d <- dual_num(3, 2)
  result <- -d
  expect_equal(primal(result), -3)
  expect_equal(tangent(result), -2)
})

test_that("dual multiplication works", {
  # d/dx (x * y) = y (when differentiating w.r.t. x)
  d1 <- dual_num(2, 1)  # x = 2, dx = 1
  d2 <- dual_num(3, 0)  # y = 3, dy = 0
  result <- d1 * d2
  expect_equal(primal(result), 6)
  expect_equal(tangent(result), 3)  # d(xy)/dx = y = 3
})

test_that("dual division works", {
  d1 <- dual_num(6, 1)  # x = 6, dx = 1
  d2 <- dual_num(2, 0)  # y = 2, dy = 0
  result <- d1 / d2
  expect_equal(primal(result), 3)
  expect_equal(tangent(result), 0.5)  # d(x/y)/dx = 1/y = 0.5
})

test_that("dual power works", {
  # f(x) = x^2 at x = 3
  d <- dual_num(3, 1)
  result <- d^dual_num(2, 0)
  expect_equal(primal(result), 9)
  expect_equal(tangent(result), 6)  # d(x^2)/dx = 2x = 6
})

test_that("dual log works", {
  d <- dual_num(2, 1)
  result <- log(d)
  expect_equal(primal(result), log(2))
  expect_equal(tangent(result), 0.5)  # d(log(x))/dx = 1/x = 0.5
})

test_that("dual exp works", {
  d <- dual_num(1, 1)
  result <- exp(d)
  expect_equal(primal(result), exp(1))
  expect_equal(tangent(result), exp(1))  # d(exp(x))/dx = exp(x)
})

test_that("dual sqrt works", {
  d <- dual_num(4, 1)
  result <- sqrt(d)
  expect_equal(primal(result), 2)
  expect_equal(tangent(result), 0.25)  # d(sqrt(x))/dx = 1/(2*sqrt(x)) = 0.25
})

test_that("dual sin/cos work", {
  d <- dual_num(pi / 4, 1)

  sin_result <- sin(d)
  expect_equal(primal(sin_result), sin(pi / 4))
  expect_equal(tangent(sin_result), cos(pi / 4))

  cos_result <- cos(d)
  expect_equal(primal(cos_result), cos(pi / 4))
  expect_equal(tangent(cos_result), -sin(pi / 4))
})

test_that("dual lgamma works", {
  d <- dual_num(3, 1)
  result <- lgamma(d)
  expect_equal(primal(result), lgamma(3))
  expect_equal(tangent(result), digamma(3))
})

test_that("dual sigmoid works", {
  d <- dual_num(0, 1)
  result <- sigmoid(d)
  expect_equal(primal(result), 0.5)
  expect_equal(tangent(result), 0.25)  # sigmoid'(0) = 0.5 * 0.5 = 0.25
})

test_that("forward-mode computes correct derivative for composite function", {
  # f(x) = x^2 + 2x + 1 = (x+1)^2
  # f'(x) = 2x + 2 = 2(x+1)
  # At x = 3: f(3) = 16, f'(3) = 8
  x <- dual_num(3, 1)
  y <- x^dual_num(2, 0) + dual_num(2, 0) * x + dual_num(1, 0)
  expect_equal(primal(y), 16)
  expect_equal(tangent(y), 8)
})

test_that("forward-mode with value objects works", {
  # f(x) = x^2 at x = 3
  # Primal and tangent are both value objects
  x_primal <- val(3)
  x_tangent <- val(1)
  x <- dual$new(x_primal, x_tangent)

  two <- dual$new(val(2), val(0))
  result <- x^two

  expect_equal(data(primal(result)), 9)
  expect_equal(data(tangent(result)), 6)
})

test_that("dual operations work with vectors (issue: && with vectors)", {
  # Regression test: vector / dual then ^dual should not error
  # The issue was is.numeric(da) && da == 0 fails when da is a vector

  # Test vector divided by dual, then raised to dual power
  t_vec <- c(1, 2, 3, 4, 5)
  sigma <- dual_num(2, 1)  # sigma with tangent
  k <- dual_num(3, 0)      # k without tangent

  # This should not error with "length = 5 in coercion to logical(1)"
  result <- (t_vec / sigma)^k

  # Check values: (t/2)^3
  expect_equal(primal(result), (t_vec / 2)^3)

  # Check tangent: d/d_sigma [(t/sigma)^k] = -k * t * (t/sigma)^(k-1) / sigma^2
  # = -3 * t * (t/2)^2 / 4 = -3 * t^3 / 16
  expected_tangent <- -3 * t_vec * (t_vec / 2)^2 / 4
  expect_equal(tangent(result), expected_tangent)
})

test_that("dual arithmetic with vector primals works", {
  # Test that basic operations work when primal is a vector
  vec <- c(1, 2, 3)
  d <- dual_num(vec, 0)

  # Addition
  result <- d + dual_num(1, 0)
  expect_equal(primal(result), c(2, 3, 4))

  # Multiplication
  result <- d * dual_num(2, 0)
  expect_equal(primal(result), c(2, 4, 6))

  # Power
  result <- d^dual_num(2, 0)
  expect_equal(primal(result), c(1, 4, 9))

  # Log
  result <- log(d)
  expect_equal(primal(result), log(vec))
})
