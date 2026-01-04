# Tests verifying gradients and Hessians against known analytical solutions
# These tests use functions with closed-form derivatives to validate correctness

# =============================================================================
# GRADIENT TESTS - Known analytical gradients
# =============================================================================

test_that("gradient of x^n matches analytical n*x^(n-1)", {
  # f(x) = x^n, f'(x) = n*x^(n-1)
  test_cases <- list(
    list(x = 2, n = 2, expected = 4),    # 2*2^1 = 4
    list(x = 3, n = 3, expected = 27),   # 3*3^2 = 27
    list(x = 2, n = 4, expected = 32),   # 4*2^3 = 32
    list(x = 5, n = 1, expected = 1),    # 1*5^0 = 1
    list(x = 4, n = 0.5, expected = 0.25) # 0.5*4^(-0.5) = 0.25
  )

  for (tc in test_cases) {
    x <- val(tc$x)
    y <- x^tc$n
    backward(y)
    expect_equal(grad(x), tc$expected, tolerance = 1e-10,
                 label = sprintf("d/dx(x^%g) at x=%g", tc$n, tc$x))
  }
})

test_that("gradient of exp(ax) matches analytical a*exp(ax)", {
  # f(x) = exp(ax), f'(x) = a*exp(ax)
  test_cases <- list(
    list(x = 0, a = 1, expected = 1),      # 1*exp(0) = 1
    list(x = 1, a = 2, expected = 2*exp(2)), # 2*exp(2)
    list(x = 2, a = -1, expected = -exp(-2)) # -1*exp(-2)
  )

  for (tc in test_cases) {
    x <- val(tc$x)
    y <- exp(tc$a * x)
    backward(y)
    expect_equal(grad(x), tc$expected, tolerance = 1e-10,
                 label = sprintf("d/dx(exp(%g*x)) at x=%g", tc$a, tc$x))
  }
})

test_that("gradient of log(x) matches analytical 1/x", {
  # f(x) = log(x), f'(x) = 1/x
  test_cases <- c(0.5, 1, 2, 5, 10)

  for (x_val in test_cases) {
    x <- val(x_val)
    y <- log(x)
    backward(y)
    expect_equal(grad(x), 1/x_val, tolerance = 1e-10,
                 label = sprintf("d/dx(log(x)) at x=%g", x_val))
  }
})

test_that("gradient of sin(x) matches analytical cos(x)", {
  # f(x) = sin(x), f'(x) = cos(x)
  test_cases <- c(0, pi/6, pi/4, pi/3, pi/2, pi)

  for (x_val in test_cases) {
    x <- val(x_val)
    y <- sin(x)
    backward(y)
    expect_equal(grad(x), cos(x_val), tolerance = 1e-10,
                 label = sprintf("d/dx(sin(x)) at x=%g", x_val))
  }
})

test_that("gradient of cos(x) matches analytical -sin(x)", {
  # f(x) = cos(x), f'(x) = -sin(x)
  test_cases <- c(0, pi/6, pi/4, pi/3, pi/2, pi)

  for (x_val in test_cases) {
    x <- val(x_val)
    y <- cos(x)
    backward(y)
    expect_equal(grad(x), -sin(x_val), tolerance = 1e-10,
                 label = sprintf("d/dx(cos(x)) at x=%g", x_val))
  }
})

test_that("gradient of tanh(x) matches analytical sech^2(x) = 1-tanh^2(x)", {
  # f(x) = tanh(x), f'(x) = 1 - tanh^2(x)
  test_cases <- c(-2, -1, 0, 1, 2)

  for (x_val in test_cases) {
    x <- val(x_val)
    y <- tanh(x)
    backward(y)
    expected <- 1 - tanh(x_val)^2
    expect_equal(grad(x), expected, tolerance = 1e-10,
                 label = sprintf("d/dx(tanh(x)) at x=%g", x_val))
  }
})

test_that("gradient of sigmoid(x) matches analytical sigmoid(x)*(1-sigmoid(x))", {
  # f(x) = 1/(1+exp(-x)), f'(x) = f(x)*(1-f(x))
  test_cases <- c(-3, -1, 0, 1, 3)

  for (x_val in test_cases) {
    x <- val(x_val)
    y <- sigmoid(x)
    backward(y)
    s <- 1 / (1 + exp(-x_val))
    expected <- s * (1 - s)
    expect_equal(grad(x), expected, tolerance = 1e-10,
                 label = sprintf("d/dx(sigmoid(x)) at x=%g", x_val))
  }
})

test_that("gradient of sqrt(x) matches analytical 1/(2*sqrt(x))", {
  # f(x) = sqrt(x), f'(x) = 1/(2*sqrt(x))
  test_cases <- c(0.25, 1, 4, 9, 16)

  for (x_val in test_cases) {
    x <- val(x_val)
    y <- sqrt(x)
    backward(y)
    expect_equal(grad(x), 1/(2*sqrt(x_val)), tolerance = 1e-10,
                 label = sprintf("d/dx(sqrt(x)) at x=%g", x_val))
  }
})

test_that("gradient of lgamma(x) matches analytical digamma(x)", {
  # f(x) = lgamma(x), f'(x) = digamma(x)
  test_cases <- c(1, 2, 3, 5, 10)

  for (x_val in test_cases) {
    x <- val(x_val)
    y <- lgamma(x)
    backward(y)
    expect_equal(grad(x), digamma(x_val), tolerance = 1e-10,
                 label = sprintf("d/dx(lgamma(x)) at x=%g", x_val))
  }
})

test_that("gradient of digamma(x) matches analytical trigamma(x)", {
  # f(x) = digamma(x), f'(x) = trigamma(x)
  test_cases <- c(1, 2, 3, 5, 10)

  for (x_val in test_cases) {
    x <- val(x_val)
    y <- digamma(x)
    backward(y)
    expect_equal(grad(x), trigamma(x_val), tolerance = 1e-10,
                 label = sprintf("d/dx(digamma(x)) at x=%g", x_val))
  }
})

test_that("gradient of log1p(x) matches analytical 1/(1+x)", {
  # f(x) = log(1+x), f'(x) = 1/(1+x)
  test_cases <- c(0, 0.5, 1, 2, 5)

  for (x_val in test_cases) {
    x <- val(x_val)
    y <- log1p(x)
    backward(y)
    expect_equal(grad(x), 1/(1+x_val), tolerance = 1e-10,
                 label = sprintf("d/dx(log1p(x)) at x=%g", x_val))
  }
})

test_that("chain rule: gradient of f(g(x)) matches f'(g(x))*g'(x)", {
  # f(x) = exp(x^2), f'(x) = 2x*exp(x^2)
  test_cases <- c(-2, -1, 0, 1, 2)

  for (x_val in test_cases) {
    x <- val(x_val)
    y <- exp(x^2)
    backward(y)
    expected <- 2 * x_val * exp(x_val^2)
    expect_equal(grad(x), expected, tolerance = 1e-10,
                 label = sprintf("d/dx(exp(x^2)) at x=%g", x_val))
  }
})

test_that("product rule: gradient of f(x)*g(x) matches f'g + fg'", {
  # f(x) = x*exp(x), f'(x) = exp(x) + x*exp(x) = (1+x)*exp(x)
  test_cases <- c(-2, -1, 0, 1, 2)

  for (x_val in test_cases) {
    x <- val(x_val)
    y <- x * exp(x)
    backward(y)
    expected <- (1 + x_val) * exp(x_val)
    expect_equal(grad(x), expected, tolerance = 1e-10,
                 label = sprintf("d/dx(x*exp(x)) at x=%g", x_val))
  }
})

test_that("quotient rule: gradient of f(x)/g(x) matches (f'g - fg')/g^2", {
  # f(x) = x/(1+x^2), f'(x) = (1*(1+x^2) - x*2x)/(1+x^2)^2 = (1-x^2)/(1+x^2)^2
  test_cases <- c(-2, -1, 0, 1, 2)

  for (x_val in test_cases) {
    x <- val(x_val)
    y <- x / (1 + x^2)
    backward(y)
    expected <- (1 - x_val^2) / (1 + x_val^2)^2
    expect_equal(grad(x), expected, tolerance = 1e-10,
                 label = sprintf("d/dx(x/(1+x^2)) at x=%g", x_val))
  }
})

test_that("multivariate gradient of f(x,y) = x^2 + y^2 matches (2x, 2y)", {
  test_cases <- list(
    list(x = 1, y = 2),
    list(x = 3, y = 4),
    list(x = -1, y = 5)
  )

  for (tc in test_cases) {
    x <- val(tc$x)
    y <- val(tc$y)
    z <- x^2 + y^2
    backward(z)
    expect_equal(grad(x), 2*tc$x, tolerance = 1e-10)
    expect_equal(grad(y), 2*tc$y, tolerance = 1e-10)
  }
})

test_that("multivariate gradient of f(x,y) = x*y matches (y, x)", {
  test_cases <- list(
    list(x = 2, y = 3),
    list(x = 5, y = 7),
    list(x = -2, y = 4)
  )

  for (tc in test_cases) {
    x <- val(tc$x)
    y <- val(tc$y)
    z <- x * y
    backward(z)
    expect_equal(grad(x), tc$y, tolerance = 1e-10)
    expect_equal(grad(y), tc$x, tolerance = 1e-10)
  }
})

test_that("multivariate gradient of Rosenbrock function is correct", {
  # f(x,y) = (1-x)^2 + 100*(y-x^2)^2
  # df/dx = -2*(1-x) + 100*2*(y-x^2)*(-2x) = -2*(1-x) - 400*x*(y-x^2)
  # df/dy = 100*2*(y-x^2) = 200*(y-x^2)

  test_cases <- list(
    list(x = 0, y = 0),
    list(x = 1, y = 1),  # minimum
    list(x = 2, y = 4),
    list(x = -1, y = 1)
  )

  for (tc in test_cases) {
    x <- val(tc$x)
    y <- val(tc$y)
    z <- (1 - x)^2 + 100 * (y - x^2)^2
    backward(z)

    expected_dx <- -2*(1-tc$x) - 400*tc$x*(tc$y - tc$x^2)
    expected_dy <- 200*(tc$y - tc$x^2)

    expect_equal(grad(x), expected_dx, tolerance = 1e-9,
                 label = sprintf("df/dx at (%g,%g)", tc$x, tc$y))
    expect_equal(grad(y), expected_dy, tolerance = 1e-9,
                 label = sprintf("df/dy at (%g,%g)", tc$x, tc$y))
  }
})

# =============================================================================
# HESSIAN TESTS - Known analytical Hessians
# =============================================================================

test_that("Hessian of f(x) = x^2 is [[2]]", {
  # f(x) = x^2, f'(x) = 2x, f''(x) = 2
  f <- function(p) p[[1]]^2

  for (x_val in c(-2, 0, 1, 3)) {
    H <- hessian(f, list(val(x_val)))
    expect_equal(H[1,1], 2, tolerance = 1e-10,
                 label = sprintf("d2/dx2(x^2) at x=%g", x_val))
  }
})

test_that("Hessian of f(x) = x^3 is [[6x]]", {
  # f(x) = x^3, f'(x) = 3x^2, f''(x) = 6x
  f <- function(p) p[[1]]^3

  for (x_val in c(-2, 1, 3)) {
    H <- hessian(f, list(val(x_val)))
    expect_equal(H[1,1], 6*x_val, tolerance = 1e-10,
                 label = sprintf("d2/dx2(x^3) at x=%g", x_val))
  }
})

test_that("Hessian of f(x) = exp(x) is [[exp(x)]]", {
  # f(x) = exp(x), f'(x) = exp(x), f''(x) = exp(x)
  f <- function(p) exp(p[[1]])

  for (x_val in c(-1, 0, 1, 2)) {
    H <- hessian(f, list(val(x_val)))
    expect_equal(H[1,1], exp(x_val), tolerance = 1e-10,
                 label = sprintf("d2/dx2(exp(x)) at x=%g", x_val))
  }
})

test_that("Hessian of f(x) = log(x) is [[-1/x^2]]", {
  # f(x) = log(x), f'(x) = 1/x, f''(x) = -1/x^2
  f <- function(p) log(p[[1]])

  for (x_val in c(0.5, 1, 2, 5)) {
    H <- hessian(f, list(val(x_val)))
    expect_equal(H[1,1], -1/x_val^2, tolerance = 1e-10,
                 label = sprintf("d2/dx2(log(x)) at x=%g", x_val))
  }
})

test_that("Hessian of f(x) = sin(x) is [[-sin(x)]]", {
  # f(x) = sin(x), f'(x) = cos(x), f''(x) = -sin(x)
  f <- function(p) sin(p[[1]])

  for (x_val in c(0, pi/4, pi/2, pi)) {
    H <- hessian(f, list(val(x_val)))
    expect_equal(H[1,1], -sin(x_val), tolerance = 1e-10,
                 label = sprintf("d2/dx2(sin(x)) at x=%g", x_val))
  }
})

test_that("Hessian of f(x) = cos(x) is [[-cos(x)]]", {
  # f(x) = cos(x), f'(x) = -sin(x), f''(x) = -cos(x)
  f <- function(p) cos(p[[1]])

  for (x_val in c(0, pi/4, pi/2, pi)) {
    H <- hessian(f, list(val(x_val)))
    expect_equal(H[1,1], -cos(x_val), tolerance = 1e-10,
                 label = sprintf("d2/dx2(cos(x)) at x=%g", x_val))
  }
})

test_that("Hessian of f(x,y) = x^2 + y^2 is [[2,0],[0,2]]", {
  # f(x,y) = x^2 + y^2
  # df/dx = 2x, df/dy = 2y
  # d2f/dx2 = 2, d2f/dxdy = 0, d2f/dy2 = 2
  f <- function(p) p[[1]]^2 + p[[2]]^2

  test_cases <- list(
    list(x = 1, y = 2),
    list(x = 0, y = 0),
    list(x = -3, y = 4)
  )

  for (tc in test_cases) {
    H <- hessian(f, list(val(tc$x), val(tc$y)))
    expect_equal(H[1,1], 2, tolerance = 1e-10)
    expect_equal(H[1,2], 0, tolerance = 1e-10)
    expect_equal(H[2,1], 0, tolerance = 1e-10)
    expect_equal(H[2,2], 2, tolerance = 1e-10)
  }
})

test_that("Hessian of f(x,y) = x*y is [[0,1],[1,0]]", {
  # f(x,y) = x*y
  # df/dx = y, df/dy = x
  # d2f/dx2 = 0, d2f/dxdy = 1, d2f/dy2 = 0
  f <- function(p) p[[1]] * p[[2]]

  test_cases <- list(
    list(x = 2, y = 3),
    list(x = 0, y = 5),
    list(x = -1, y = 4)
  )

  for (tc in test_cases) {
    H <- hessian(f, list(val(tc$x), val(tc$y)))
    expect_equal(H[1,1], 0, tolerance = 1e-10)
    expect_equal(H[1,2], 1, tolerance = 1e-10)
    expect_equal(H[2,1], 1, tolerance = 1e-10)
    expect_equal(H[2,2], 0, tolerance = 1e-10)
  }
})

test_that("Hessian of f(x,y) = x^2*y + x*y^2 is correct", {
  # f(x,y) = x^2*y + x*y^2
  # df/dx = 2xy + y^2
  # df/dy = x^2 + 2xy
  # d2f/dx2 = 2y
  # d2f/dxdy = 2x + 2y
  # d2f/dy2 = 2x
  f <- function(p) {
    x <- p[[1]]
    y <- p[[2]]
    x^2 * y + x * y^2
  }

  test_cases <- list(
    list(x = 1, y = 2),
    list(x = 3, y = 1),
    list(x = 2, y = 2)
  )

  for (tc in test_cases) {
    H <- hessian(f, list(val(tc$x), val(tc$y)))
    expect_equal(H[1,1], 2*tc$y, tolerance = 1e-10,
                 label = sprintf("d2f/dx2 at (%g,%g)", tc$x, tc$y))
    expect_equal(H[1,2], 2*tc$x + 2*tc$y, tolerance = 1e-10,
                 label = sprintf("d2f/dxdy at (%g,%g)", tc$x, tc$y))
    expect_equal(H[2,1], 2*tc$x + 2*tc$y, tolerance = 1e-10,
                 label = sprintf("d2f/dydx at (%g,%g)", tc$x, tc$y))
    expect_equal(H[2,2], 2*tc$x, tolerance = 1e-10,
                 label = sprintf("d2f/dy2 at (%g,%g)", tc$x, tc$y))
  }
})

test_that("Hessian of quadratic form f(x) = x'Ax is 2A for symmetric A", {
  # f(x,y) = ax^2 + 2bxy + cy^2 = [x,y] * [[a,b],[b,c]] * [x,y]'
  # Hessian = [[2a, 2b], [2b, 2c]]

  # Test case: f(x,y) = 3x^2 + 4xy + 5y^2
  # A = [[3,2],[2,5]], Hessian = [[6,4],[4,10]]
  f <- function(p) {
    x <- p[[1]]
    y <- p[[2]]
    3*x^2 + 4*x*y + 5*y^2
  }

  H <- hessian(f, list(val(1), val(2)))
  expect_equal(H[1,1], 6, tolerance = 1e-10)
  expect_equal(H[1,2], 4, tolerance = 1e-10)
  expect_equal(H[2,1], 4, tolerance = 1e-10)
  expect_equal(H[2,2], 10, tolerance = 1e-10)
})

test_that("Hessian of Rosenbrock function is correct", {
  # f(x,y) = (1-x)^2 + 100*(y-x^2)^2
  # df/dx = -2*(1-x) - 400*x*(y-x^2)
  # df/dy = 200*(y-x^2)
  # d2f/dx2 = 2 - 400*(y-x^2) + 800*x^2 = 2 - 400*y + 1200*x^2
  # d2f/dxdy = -400*x
  # d2f/dy2 = 200
  f <- function(p) {
    x <- p[[1]]
    y <- p[[2]]
    (1 - x)^2 + 100 * (y - x^2)^2
  }

  test_cases <- list(
    list(x = 0, y = 0),
    list(x = 1, y = 1),  # minimum
    list(x = 2, y = 4)
  )

  for (tc in test_cases) {
    H <- hessian(f, list(val(tc$x), val(tc$y)))

    expected_xx <- 2 - 400*tc$y + 1200*tc$x^2
    expected_xy <- -400*tc$x
    expected_yy <- 200

    expect_equal(H[1,1], expected_xx, tolerance = 1e-8,
                 label = sprintf("d2f/dx2 at (%g,%g)", tc$x, tc$y))
    expect_equal(H[1,2], expected_xy, tolerance = 1e-8,
                 label = sprintf("d2f/dxdy at (%g,%g)", tc$x, tc$y))
    expect_equal(H[2,1], expected_xy, tolerance = 1e-8,
                 label = sprintf("d2f/dydx at (%g,%g)", tc$x, tc$y))
    expect_equal(H[2,2], expected_yy, tolerance = 1e-8,
                 label = sprintf("d2f/dy2 at (%g,%g)", tc$x, tc$y))
  }
})

test_that("Hessian of f(x,y) = exp(x+y) is exp(x+y)*[[1,1],[1,1]]", {
  # f(x,y) = exp(x+y)
  # df/dx = exp(x+y), df/dy = exp(x+y)
  # d2f/dx2 = exp(x+y), d2f/dxdy = exp(x+y), d2f/dy2 = exp(x+y)
  f <- function(p) exp(p[[1]] + p[[2]])

  test_cases <- list(
    list(x = 0, y = 0),
    list(x = 1, y = -1),
    list(x = 0.5, y = 0.5)
  )

  for (tc in test_cases) {
    H <- hessian(f, list(val(tc$x), val(tc$y)))
    e <- exp(tc$x + tc$y)
    expect_equal(H[1,1], e, tolerance = 1e-10)
    expect_equal(H[1,2], e, tolerance = 1e-10)
    expect_equal(H[2,1], e, tolerance = 1e-10)
    expect_equal(H[2,2], e, tolerance = 1e-10)
  }
})

test_that("Hessian of f(x,y) = log(x*y) is [[-1/x^2, 0], [0, -1/y^2]]", {
  # f(x,y) = log(x*y) = log(x) + log(y)
  # df/dx = 1/x, df/dy = 1/y
  # d2f/dx2 = -1/x^2, d2f/dxdy = 0, d2f/dy2 = -1/y^2
  f <- function(p) log(p[[1]] * p[[2]])

  test_cases <- list(
    list(x = 1, y = 1),
    list(x = 2, y = 3),
    list(x = 0.5, y = 4)
  )

  for (tc in test_cases) {
    H <- hessian(f, list(val(tc$x), val(tc$y)))
    expect_equal(H[1,1], -1/tc$x^2, tolerance = 1e-10)
    expect_equal(H[1,2], 0, tolerance = 1e-10)
    expect_equal(H[2,1], 0, tolerance = 1e-10)
    expect_equal(H[2,2], -1/tc$y^2, tolerance = 1e-10)
  }
})

test_that("Hessian symmetry holds for various functions", {
  # Hessian should always be symmetric (Schwarz's theorem)
  functions <- list(
    function(p) p[[1]]^3 * p[[2]]^2,
    function(p) sin(p[[1]]) * cos(p[[2]]),
    function(p) exp(p[[1]] * p[[2]]),
    function(p) log(p[[1]]^2 + p[[2]]^2 + 1)
  )

  for (f in functions) {
    H <- hessian(f, list(val(1.5), val(2.5)))
    expect_equal(H[1,2], H[2,1], tolerance = 1e-9,
                 label = "Hessian symmetry")
  }
})

test_that("3x3 Hessian of f(x,y,z) = x^2 + y^2 + z^2 + xy + yz is correct", {
  # f = x^2 + y^2 + z^2 + xy + yz
  # df/dx = 2x + y
  # df/dy = 2y + x + z
  # df/dz = 2z + y
  # Hessian = [[2, 1, 0], [1, 2, 1], [0, 1, 2]]
  f <- function(p) {
    x <- p[[1]]
    y <- p[[2]]
    z <- p[[3]]
    x^2 + y^2 + z^2 + x*y + y*z
  }

  H <- hessian(f, list(val(1), val(2), val(3)))

  expected <- matrix(c(2, 1, 0,
                       1, 2, 1,
                       0, 1, 2), nrow = 3, byrow = TRUE)

  expect_equal(H, expected, tolerance = 1e-10)
})

# =============================================================================
# STATISTICAL FUNCTION TESTS
# =============================================================================

test_that("gradient of normal log-likelihood is score function", {
  # loglik = -n/2*log(2*pi) - n*log(sigma) - sum((x-mu)^2)/(2*sigma^2)
  # d/dmu = sum(x-mu)/sigma^2
  # d/dsigma = -n/sigma + sum((x-mu)^2)/sigma^3

  set.seed(42)
  x <- rnorm(50, mean = 5, sd = 2)
  n <- length(x)

  mu <- val(5)
  sigma <- val(2)
  ll <- loglik_normal(mu, sigma, x)
  backward(ll)

  # Analytical score
  expected_dmu <- sum(x - 5) / 4
  expected_dsigma <- -n/2 + sum((x - 5)^2) / 8

  expect_equal(grad(mu), expected_dmu, tolerance = 1e-10)
  expect_equal(grad(sigma), expected_dsigma, tolerance = 1e-10)
})

test_that("gradient of exponential log-likelihood is score function", {
  # loglik = n*log(lambda) - lambda*sum(x)
  # d/dlambda = n/lambda - sum(x)

  set.seed(42)
  x <- rexp(50, rate = 2)
  n <- length(x)

  lambda <- val(2)
  ll <- loglik_exponential(lambda, x)
  backward(ll)

  expected <- n/2 - sum(x)
  expect_equal(grad(lambda), expected, tolerance = 1e-10)
})

test_that("gradient of Poisson log-likelihood is score function", {
  # loglik = sum(x)*log(lambda) - n*lambda - sum(lgamma(x+1))
  # d/dlambda = sum(x)/lambda - n

  set.seed(42)
  x <- rpois(50, lambda = 3)
  n <- length(x)

  lambda <- val(3)
  ll <- loglik_poisson(lambda, x)
  backward(ll)

  expected <- sum(x)/3 - n
  expect_equal(grad(lambda), expected, tolerance = 1e-10)
})

test_that("Hessian of normal log-likelihood gives Fisher information", {
  # For normal: I(mu, sigma) = [[n/sigma^2, 0], [0, 2n/sigma^2]]
  # Hessian of loglik = -I (negative Fisher information)

  set.seed(42)
  x <- rnorm(100, mean = 5, sd = 2)
  n <- length(x)

  loglik <- function(p) loglik_normal(p[[1]], p[[2]], x)

  # Evaluate at true parameters
  H <- hessian(loglik, list(val(5), val(2)))

  # Expected negative Fisher information
  # d2/dmu2 = -n/sigma^2 = -100/4 = -25
  # d2/dmudsigma = 0 (at MLE)
  # d2/dsigma2 = -2n/sigma^2 = -200/4 = -50 (at MLE with mean-centered data)

  expect_equal(H[1,1], -n/4, tolerance = 1e-10)
  # Cross terms depend on data, but should be small at true params
})

test_that("Hessian of quadratic loss is constant", {
  # L = sum((y - (a + b*x))^2)
  # dL/da = -2*sum(y - a - b*x)
  # dL/db = -2*sum(x*(y - a - b*x))
  # d2L/da2 = 2*n
  # d2L/dadb = 2*sum(x)
  # d2L/db2 = 2*sum(x^2)

  set.seed(42)
  x <- 1:10
  y <- 2 + 3*x + rnorm(10, sd = 0.1)

  loss <- function(p) {
    a <- p[[1]]
    b <- p[[2]]
    residuals <- y - (a + b*x)
    sum(residuals^2)
  }

  H <- hessian(loss, list(val(2), val(3)))

  n <- length(x)
  expect_equal(H[1,1], 2*n, tolerance = 1e-10)
  expect_equal(H[1,2], 2*sum(x), tolerance = 1e-10)
  expect_equal(H[2,1], 2*sum(x), tolerance = 1e-10)
  expect_equal(H[2,2], 2*sum(x^2), tolerance = 1e-10)
})
