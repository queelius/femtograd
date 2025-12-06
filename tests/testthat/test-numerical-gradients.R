test_that("gradients match numerical gradients for simple functions", {
  # Helper function to compute numerical gradient
  numerical_grad <- function(f, x, h = 1e-5) {
    (f(x + h) - f(x - h)) / (2 * h)
  }

  # Test f(x) = x^2
  x_val <- 3
  x <- val(x_val)
  y <- x^2
  backward(y)

  numerical <- numerical_grad(function(t) t^2, x_val)
  expect_equal(grad(x), numerical, tolerance = 1e-4)

  # Test f(x) = log(x)
  x_val <- 2
  x <- val(x_val)
  y <- log(x)
  backward(y)

  numerical <- numerical_grad(function(t) log(t), x_val)
  expect_equal(grad(x), numerical, tolerance = 1e-4)

  # Test f(x) = exp(x)
  x_val <- 1.5
  x <- val(x_val)
  y <- exp(x)
  backward(y)

  numerical <- numerical_grad(function(t) exp(t), x_val)
  expect_equal(grad(x), numerical, tolerance = 1e-4)
})

test_that("gradients match numerical gradients for abs", {
  numerical_grad <- function(f, x, h = 1e-5) {
    (f(x + h) - f(x - h)) / (2 * h)
  }

  # Test positive value
  x_val <- 3
  x <- val(x_val)
  y <- abs(x)
  backward(y)

  numerical <- numerical_grad(function(t) abs(t), x_val)
  expect_equal(grad(x), numerical, tolerance = 1e-4)

  # Test negative value
  x_val <- -3
  x <- val(x_val)
  y <- abs(x)
  backward(y)

  numerical <- numerical_grad(function(t) abs(t), x_val)
  expect_equal(grad(x), numerical, tolerance = 1e-4)
})

test_that("gradients match numerical gradients for composite functions", {
  numerical_grad_2d <- function(f, x_val, y_val, h = 1e-5) {
    df_dx <- (f(x_val + h, y_val) - f(x_val - h, y_val)) / (2 * h)
    df_dy <- (f(x_val, y_val + h) - f(x_val, y_val - h)) / (2 * h)
    c(df_dx, df_dy)
  }

  # Test f(x, y) = x*y + x^2
  x_val <- 3
  y_val <- 4
  x <- val(x_val)
  y <- val(y_val)
  z <- x * y + x^2
  backward(z)

  numerical <- numerical_grad_2d(function(a, b) a * b + a^2, x_val, y_val)
  expect_equal(grad(x), numerical[1], tolerance = 1e-4)
  expect_equal(grad(y), numerical[2], tolerance = 1e-4)
})

test_that("L1 loss gradient is correct", {
  # L1 loss: sum(abs(predictions - targets))
  predictions <- list(val(2), val(3), val(5))
  targets <- c(1.5, 4, 4.5)

  losses <- mapply(function(pred, target) abs(pred - target),
                   predictions, targets, SIMPLIFY = FALSE)

  total_loss <- sum(losses[[1]], losses[[2]], losses[[3]])
  backward(total_loss)

  # Check gradients
  # pred[1] = 2, target[1] = 1.5, diff = 0.5 > 0, so grad = 1
  expect_equal(grad(predictions[[1]]), 1)

  # pred[2] = 3, target[2] = 4, diff = -1 < 0, so grad = -1
  expect_equal(grad(predictions[[2]]), -1)

  # pred[3] = 5, target[3] = 4.5, diff = 0.5 > 0, so grad = 1
  expect_equal(grad(predictions[[3]]), 1)
})
