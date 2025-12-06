test_that("hessian works for quadratic function", {
  # f(x, y) = x^2 + xy + y^2
  # H = [[2, 1], [1, 2]]
  loss_fn <- function(p) {
    x <- p[[1]]
    y <- p[[2]]
    x^2 + x * y + y^2
  }

  params <- list(val(1), val(2))
  H <- hessian(loss_fn, params)

  expect_equal(H[1, 1], 2, tolerance = 1e-6)
  expect_equal(H[1, 2], 1, tolerance = 1e-6)
  expect_equal(H[2, 1], 1, tolerance = 1e-6)
  expect_equal(H[2, 2], 2, tolerance = 1e-6)
})

test_that("hessian is symmetric", {
  # f(x, y, z) = x*y + y*z + x*z
  loss_fn <- function(p) {
    x <- p[[1]]
    y <- p[[2]]
    z <- p[[3]]
    x * y + y * z + x * z
  }

  params <- list(val(1), val(2), val(3))
  H <- hessian(loss_fn, params)

  # Check symmetry
  expect_equal(H[1, 2], H[2, 1], tolerance = 1e-6)
  expect_equal(H[1, 3], H[3, 1], tolerance = 1e-6)
  expect_equal(H[2, 3], H[3, 2], tolerance = 1e-6)

  # For this function, diagonal is 0, off-diagonal is 1
  expect_equal(H[1, 1], 0, tolerance = 1e-6)
  expect_equal(H[1, 2], 1, tolerance = 1e-6)
})

test_that("hessian works for exponential function", {
  # f(x) = exp(x), f''(x) = exp(x)
  loss_fn <- function(p) exp(p[[1]])

  params <- list(val(2))
  H <- hessian(loss_fn, params)

  expect_equal(H[1, 1], exp(2), tolerance = 1e-5)
})

test_that("hessian works for log function", {
  # f(x) = log(x), f''(x) = -1/x^2
  loss_fn <- function(p) log(p[[1]])

  params <- list(val(3))
  H <- hessian(loss_fn, params)

  expect_equal(H[1, 1], -1 / 9, tolerance = 1e-5)
})

test_that("gradient function returns correct values", {
  # f(x, y) = x^2 + y*2
  # grad = (2x, 2)
  loss_fn <- function(p) {
    p[[1]]^2 + p[[2]] * 2
  }

  params <- list(val(3), val(5))
  g <- gradient(loss_fn, params)

  expect_equal(g[1], 6)  # 2*3 = 6
  expect_equal(g[2], 2)
})

test_that("fisher_information is negative hessian", {
  loss_fn <- function(p) {
    x <- p[[1]]
    -x^2  # negative quadratic, hessian = -2
  }

  params <- list(val(1))
  H <- hessian(loss_fn, params)
  I <- fisher_information(loss_fn, params)

  expect_equal(I[1, 1], -H[1, 1])
})

test_that("std_errors computes correctly", {
  # For a quadratic -x^2/2, Hessian = -1
  # Variance = -1/H = 1, SE = 1
  H <- matrix(-1, 1, 1)
  se <- std_errors(H, is_loglik = TRUE)
  expect_equal(se[1], 1)
})

test_that("vcov_matrix computes correctly", {
  H <- matrix(c(-2, 0, 0, -4), 2, 2)
  vcov <- vcov_matrix(H, is_loglik = TRUE)

  # -H^-1 = [[0.5, 0], [0, 0.25]]
  expect_equal(vcov[1, 1], 0.5)
  expect_equal(vcov[2, 2], 0.25)
})
