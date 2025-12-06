test_that("loglik_exponential computes correct value", {
  x <- c(1, 2, 3)
  rate <- val(2)

  ll <- loglik_exponential(rate, x)
  # n*log(rate) - rate*sum(x) = 3*log(2) - 2*6 = 3*log(2) - 12
  expected <- 3 * log(2) - 12
  expect_equal(data(ll), expected)
})

test_that("loglik_exponential gradient is score function", {
  x <- c(1, 2, 3)
  rate <- val(2)

  ll <- loglik_exponential(rate, x)
  backward(ll)

  # d/d(rate) [n*log(rate) - rate*sum(x)] = n/rate - sum(x)
  expected_grad <- 3 / 2 - 6
  expect_equal(grad(rate), expected_grad)
})

test_that("loglik_normal computes correct value", {
  x <- c(0, 1, 2)
  mu <- val(1)
  sigma <- val(1)

  ll <- loglik_normal(mu, sigma, x)

  # Manual calculation
  n <- 3
  sum_sq <- sum((x - 1)^2)  # 1 + 0 + 1 = 2
  expected <- -n / 2 * log(2 * pi) - n * log(1) - sum_sq / 2
  expect_equal(data(ll), expected, tolerance = 1e-10)
})

test_that("loglik_normal gradient w.r.t. mu is correct", {
  x <- c(0, 2)  # mean = 1
  mu <- val(1)
  sigma <- val(1)

  ll <- loglik_normal(mu, sigma, x)
  backward(ll)

  # At MLE (mu = mean(x)), score should be 0
  expect_equal(grad(mu), 0, tolerance = 1e-10)
})

test_that("loglik_poisson computes correct value", {
  x <- c(1, 2, 3)
  lambda <- val(2)

  ll <- loglik_poisson(lambda, x)

  # sum(x)*log(lambda) - n*lambda - sum(lgamma(x+1))
  expected <- sum(x) * log(2) - 3 * 2 - sum(lgamma(x + 1))
  expect_equal(data(ll), expected, tolerance = 1e-10)
})

test_that("loglik_poisson gradient is correct", {
  x <- c(1, 2, 3)
  lambda <- val(2)

  ll <- loglik_poisson(lambda, x)
  backward(ll)

  # d/d(lambda) = sum(x)/lambda - n
  expected <- sum(x) / 2 - 3
  expect_equal(grad(lambda), expected)
})

test_that("loglik_bernoulli computes correct value", {
  x <- c(1, 0, 1, 1)  # 3 successes, 1 failure
  p <- val(0.5)

  ll <- loglik_bernoulli(p, x)

  # 3*log(0.5) + 1*log(0.5) = 4*log(0.5)
  expected <- 4 * log(0.5)
  expect_equal(data(ll), expected)
})

test_that("loglik_binomial gradient is correct", {
  x <- c(3, 4, 5)  # successes
  n <- c(10, 10, 10)  # trials
  p <- val(0.4)

  ll <- loglik_binomial(p, x, size = n)
  backward(ll)

  # d/dp = sum(x)/p - sum(n-x)/(1-p)
  expected <- sum(x) / 0.4 - sum(n - x) / 0.6
  expect_equal(grad(p), expected, tolerance = 1e-10)
})

test_that("loglik_gamma computes correct value", {
  x <- c(1, 2, 3)
  shape <- val(2)
  rate <- val(1)

  ll <- loglik_gamma(shape, rate, x)

  # n*shape*log(rate) - n*lgamma(shape) + (shape-1)*sum(log(x)) - rate*sum(x)
  n <- 3
  expected <- n * 2 * log(1) - n * lgamma(2) + (2 - 1) * sum(log(x)) - 1 * sum(x)
  expect_equal(data(ll), expected, tolerance = 1e-10)
})

test_that("loglik_beta is defined on (0,1) interval", {
  x <- c(0.2, 0.5, 0.8)
  alpha <- val(2)
  beta <- val(3)

  ll <- loglik_beta(alpha, beta, x)

  # Just check it computes without error and is finite

  expect_true(is.finite(data(ll)))
})

test_that("MLE for exponential is 1/mean(x)", {
  set.seed(42)
  x <- rexp(100, rate = 2)

  # At MLE, gradient should be zero
  mle_rate <- 1 / mean(x)
  rate <- val(mle_rate)
  ll <- loglik_exponential(rate, x)
  backward(ll)

  expect_equal(grad(rate), 0, tolerance = 1e-10)
})
