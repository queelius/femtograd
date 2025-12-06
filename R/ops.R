
#' Addition for value objects
#'
#' @param e1 A value object
#' @param e2 A value object or a scalar
#'
#' @return A new value object representing the sum
#' @export
`+.value` <- function(e1, e2)
{
  if (!is_value(e1))
    e1 <- val(e1)
  if (!is_value(e2))
    e2 <- val(e2)
  out <- value$new(e1$data + e2$data, list(e1, e2))
  out$backward_fn <- function()
  {
    e1$grad <<- e1$grad + out$grad
    e2$grad <<- e2$grad + out$grad
  }
  out
}

#' Multiplication for value objects
#'
#' @param e1 A value object
#' @param e2 A value object or a scalar
#'
#' @return A new value object representing the product
#' @export
`*.value` <- function(e1, e2)
{
  if (!is_value(e1))
    e1 <- val(e1)
  if (!is_value(e2))
    e2 <- val(e2)
  out <- value$new(e1$data * e2$data, list(e1, e2))
  out$backward_fn <- function()
  {
    e1$grad <<- e1$grad + e2$data * out$grad
    e2$grad <<- e2$grad + e1$data * out$grad
  }
  out
}

#' Power operation for value objects.
#'
#' @param b A value object as the base
#' @param e A value object, a scalar, or a numeric constant as the exponent
#'
#' @return A new value object representing the power b^e
#' @export
`^.value` <- function(b, e)
{
  # Handle scalar base (e.g., 2^val(3))
  if (!is_value(b))
    b <- val(b)

  if (is_value(e))
  {
    out <- value$new(b$data^e$data, list(b, e))
    out$backward_fn <- function()
    {
      b$grad <<- b$grad + e$data * b$data^(e$data - 1) * out$grad
      e$grad <<- e$grad + log(b$data) * b$data^e$data * out$grad
    }
    return(out)
  }
  else
  {
    out <- value$new(b$data^e, list(b))
    out$backward_fn <- function()
    {
      b$grad <<- b$grad + e * b$data^(e - 1) * out$grad
    }
    return(out)
  }
}

#' Summation for value objects
#'
#' Overrides base::sum to work with value objects. Supports mixed lists
#' of value objects and numeric scalars.
#'
#' @param ... value objects and/or numeric scalars to sum
#' @param na.rm Logical, whether to remove NA values (passed to base sum for scalars)
#'
#' @return A new value object representing the sum
#' @export
sum.value <- function(..., na.rm = FALSE)
{
  values <- list(...)
  if (length(values) == 0)
    return(val(0))

  # Extract value objects and scalars
  value_objs <- Filter(is_value, values)
  scalars <- Filter(Negate(is_value), values)

  # Compute total
  total <- sum(sapply(value_objs, function(v) v$data))
  if (length(scalars) > 0)
    total <- total + sum(unlist(scalars), na.rm = na.rm)

  out <- value$new(total, value_objs)
  out$backward_fn <- function()
  {
    for (child in out$prev)
      child$grad <- child$grad + out$grad
  }
  out
}

#' Mean for value objects
#'
#' Overrides base::mean to work with value objects. Computes the arithmetic
#' mean and properly propagates gradients (each input gets gradient / n).
#'
#' @param x A value object or vector of value objects
#' @param ... Additional arguments (for compatibility with base::mean)
#'
#' @return A new value object representing the mean
#' @export
mean.value <- function(x, ...)
{
  # Handle single value object
  if (is_value(x) && is.null(x$prev))
    return(x)

  # Handle vector/list of values
  values <- if (is.list(x)) x else list(x)
  value_objs <- Filter(is_value, values)

  if (length(value_objs) == 0)
    return(val(0))

  n <- length(value_objs)
  avg <- sum(sapply(value_objs, function(v) v$data)) / n

  out <- value$new(avg, value_objs)
  out$backward_fn <- function()
  {
    for (child in out$prev)
      child$grad <- child$grad + out$grad / n
  }
  out
}


#' Natural logarithm for value objects
#'
#' @param x A value object
#'
#' @return A new value object representing the natural logarithm of x
#' @export
log.value <- function(x)
{
  out <- value$new(log(x$data), list(x))
  out$backward_fn <- function()
  {
    x$grad <<- x$grad + (1 / x$data) * out$grad
  }
  out
}

#' Exponential function for value objects
#'
#' @param x A value object
#'
#' @return A new value object representing the exponential of x
#' @export
exp.value <- function(x)
{
  out <- value$new(exp(x$data), list(x))
  out$backward_fn <- function()
  {
    x$grad <<- x$grad + (exp(x$data)) * out$grad
  }
  out
}

#' Division for value objects
#'
#' @param x A value object (numerator)
#' @param y A value object or a scalar (denominator)
#'
#' @return A new value object representing the division x / y
#' @export
`/.value` <- function(x, y)
{
  if (!is_value(x))
    x <- val(x)
  if (!is_value(y))
    y <- val(y)
  out <- value$new(x$data / y$data, list(x, y))
  out$backward_fn <- function()
  {
    x$grad <<- x$grad + (1 / y$data) * out$grad
    y$grad <<- y$grad - (x$data / (y$data^2)) * out$grad
  }
  out
}

#' Square root for value objects
#'
#' @param x A value object
#'
#' @return A new value object representing the square root of x
#' @export
sqrt.value <- function(x)
{
  out <- value$new(sqrt(x$data), list(x))
  out$backward_fn <- function()
  {
    x$grad <<- x$grad + (1 / (2 * sqrt(x$data))) * out$grad
  }
  return(out)
}


#' Subtraction for value objects
#'
#' @param x A value object (minuend)
#' @param y A value object or a scalar (subtrahend)
#'
#' @return A new value object representing the subtraction x - y
#' @export
`-.value` <- function(x, y = NULL)
{
  # Handle unary negation
  if (missing(y)) {
    if (!is_value(x))
      x <- val(x)
    out <- value$new(-x$data, list(x))
    out$backward_fn <- function() {
      x$grad <<- x$grad - out$grad
    }
    return(out)
  }

  # Handle binary subtraction
  if (!is_value(x))
    x <- val(x)
  if (!is_value(y))
    y <- val(y)
  out <- value$new(x$data - y$data, list(x, y))
  out$backward_fn <- function()
  {
    x$grad <<- x$grad + out$grad
    y$grad <<- y$grad - out$grad
  }
  out
}

#' Hyperbolic tangent activation function for value objects
#'
#' @param x A value object
#'
#' @return A new value object representing tanh(x)
#' @export
tanh.value <- function(x)
{
  t <- tanh(x$data)
  out <- value$new(t, list(x))
  out$backward_fn <- function()
  {
    x$grad <<- x$grad + (1 - t^2) * out$grad
  }
  out
}

#' Sigmoid activation function for value objects
#'
#' @param x A value object
#'
#' @return A new value object representing sigmoid(x) = 1/(1+exp(-x))
#' @export
sigmoid <- function(x)
{
  UseMethod("sigmoid")
}

#' @export
sigmoid.value <- function(x)
{
  s <- 1 / (1 + exp(-x$data))
  out <- value$new(s, list(x))
  out$backward_fn <- function()
  {
    x$grad <<- x$grad + s * (1 - s) * out$grad
  }
  out
}

#' @export
sigmoid.default <- function(x)
{
  1 / (1 + exp(-x))
}

#' ReLU (Rectified Linear Unit) activation function for value objects
#'
#' @param x A value object
#'
#' @return A new value object representing max(0, x)
#' @export
relu <- function(x)
{
  UseMethod("relu")
}

#' @export
relu.value <- function(x)
{
  out <- value$new(max(0, x$data), list(x))
  out$backward_fn <- function()
  {
    x$grad <<- x$grad + (if (x$data > 0) out$grad else 0)
  }
  out
}

#' @export
relu.default <- function(x)
{
  max(0, x)
}

#' Absolute value for value objects
#'
#' @param x A value object
#'
#' @return A new value object representing abs(x)
#' @export
abs.value <- function(x)
{
  out <- value$new(abs(x$data), list(x))
  out$backward_fn <- function()
  {
    # Gradient of abs(x) is sign(x), undefined at 0 (we use 0 there)
    x$grad <<- x$grad + sign(x$data) * out$grad
  }
  out
}

#' Log-gamma function for value objects
#'
#' Computes log(gamma(x)), essential for many statistical distributions
#' (Poisson, binomial, gamma, beta, etc.)
#'
#' @param x A value object
#'
#' @return A new value object representing lgamma(x)
#' @export
lgamma.value <- function(x)
{
  out <- value$new(lgamma(x$data), list(x))
  out$backward_fn <- function()
  {
    # d/dx lgamma(x) = digamma(x) = psi(x)
    x$grad <<- x$grad + digamma(x$data) * out$grad
  }
  out
}

#' Digamma (psi) function for value objects
#'
#' Computes the digamma function psi(x) = d/dx log(gamma(x)).
#' This is the derivative of lgamma and appears in gradients of
#' many statistical distributions.
#'
#' @param x A value object
#'
#' @return A new value object representing digamma(x)
#' @export
digamma.value <- function(x)
{
  out <- value$new(digamma(x$data), list(x))
  out$backward_fn <- function()
  {
    # d/dx digamma(x) = trigamma(x)
    x$grad <<- x$grad + trigamma(x$data) * out$grad
  }
  out
}

#' Trigamma function for value objects
#'
#' Computes the trigamma function, the second derivative of lgamma.
#' Useful for Hessian computations involving gamma-related distributions.
#'
#' @param x A value object
#'
#' @return A new value object representing trigamma(x)
#' @export
trigamma.value <- function(x)
{
  out <- value$new(trigamma(x$data), list(x))
  out$backward_fn <- function()
  {
    # d/dx trigamma(x) = psigamma(x, 2)
    x$grad <<- x$grad + psigamma(x$data, 2) * out$grad
  }
  out
}

#' Log(1+x) for value objects
#'
#' Computes log(1+x) in a numerically stable way. Essential for
#' computing log-likelihoods with probabilities near 0 or 1.
#'
#' @param x A value object
#'
#' @return A new value object representing log(1+x)
#' @export
log1p.value <- function(x)
{
  out <- value$new(log1p(x$data), list(x))
  out$backward_fn <- function()
  {
    # d/dx log(1+x) = 1/(1+x)
    x$grad <<- x$grad + (1 / (1 + x$data)) * out$grad
  }
  out
}

#' Sine function for value objects
#'
#' @param x A value object
#'
#' @return A new value object representing sin(x)
#' @export
sin.value <- function(x)
{
  out <- value$new(sin(x$data), list(x))
  out$backward_fn <- function()
  {
    x$grad <<- x$grad + cos(x$data) * out$grad
  }
  out
}

#' Cosine function for value objects
#'
#' @param x A value object
#'
#' @return A new value object representing cos(x)
#' @export
cos.value <- function(x)
{
  out <- value$new(cos(x$data), list(x))
  out$backward_fn <- function()
  {
    x$grad <<- x$grad - sin(x$data) * out$grad
  }
  out
}

#' Logit function for value objects
#'
#' Computes log(p/(1-p)), the log-odds transformation.
#' Maps probabilities (0,1) to real line (-inf, inf).
#'
#' @param x A value object representing a probability
#'
#' @return A new value object representing logit(x)
#' @export
logit <- function(x)
{
  UseMethod("logit")
}

#' @export
logit.value <- function(x)
{
  out <- value$new(log(x$data / (1 - x$data)), list(x))
  out$backward_fn <- function()
  {
    # d/dx logit(x) = 1/(x(1-x))
    x$grad <<- x$grad + (1 / (x$data * (1 - x$data))) * out$grad
  }
  out
}

#' @export
logit.default <- function(x)
{
  log(x / (1 - x))
}

#' Softplus function for value objects
#'
#' Computes log(1 + exp(x)), a smooth approximation to max(0, x).
#' Useful for ensuring positivity of parameters (e.g., variances).
#'
#' @param x A value object
#'
#' @return A new value object representing softplus(x)
#' @export
softplus <- function(x)
{
  UseMethod("softplus")
}

#' @export
softplus.value <- function(x)
{
  # Numerically stable: for large x, softplus(x) ≈ x
  out <- value$new(log1p(exp(x$data)), list(x))
  out$backward_fn <- function()
  {
    # d/dx softplus(x) = sigmoid(x) = 1/(1+exp(-x))
    x$grad <<- x$grad + (1 / (1 + exp(-x$data))) * out$grad
  }
  out
}

#' @export
softplus.default <- function(x)
{
  log1p(exp(x))
}
