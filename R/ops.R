
#' Addition for value objects
#'
#' @param e1 A value object
#' @param e2 A value object or a scalar
#'
#' @return A new value object representing the sum
#' @export
`+.value` <- function(e1, e2)
{
  if (!is_value(e2))
    e2 <- val(e2)
  out <- value$new(e1$data + e2$data, list(e1, e2), '+')
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
  if (!is_value(e2))
    e2 <- val(e2)
  out <- value$new(e1$data * e2$data, list(e1, e2), '*')
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
  if (is_value(e))
  {
    out <- value$new(b$data^e$data, list(b, e), '^')
    out$backward_fn <- function()
    {
      b$grad <<- b$grad + e$data * b$data^(e$data - 1) * out$grad
      e$grad <<- e$grad + log(b$data) * b$data^e$data * out$grad
    }
    return(out)
  }
  else
  {
    out <- value$new(b$data^e, list(b), '^')
    out$backward_fn <- function()
    {
      b$grad <<- b$grad + e * b$data^(e - 1) * out$grad
    }
    return(out)
  }
}

#' Summation for a list of value objects
#'
#' @param values A list of objects that can be added
#'
#' @return A new value object representing the sum of all input value objects
#' @export
sum_values <- function(values = NULL)
{
  if (length(values) == 0)
    return(val(0))

  s <- sum(sapply(values, retrieve))
  out <- value$new(s, Filter(is_value, values), "sum_values")
  out$backward_fn <- function()
  {
    for (child in out$prev)
      child$grad <- child$grad + out$grad
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
  out <- value$new(log(x$data), list(x), 'log')
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
  out <- value$new(exp(x$data), list(x), 'exp')
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
  if (!is_value(y))
    y <- val(y)
  out <- value$new(x$data / y$data, list(x, y), '/')
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
  out <- value$new(sqrt(x$data), list(x), 'sqrt')
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
`-.value` <- function(x, y)
{
  y <- if (is_value(y)) y else val(y)
  out <- value$new(x$data - y$data, list(x, y), '-')
  out$backward_fn <- function()
  {
    x$grad <<- x$grad + out$grad
    y$grad <<- y$grad - out$grad
  }
  out
}
