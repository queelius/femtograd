
#' Addition for value objects
#'
#' Element-wise addition of two value objects or a value and numeric.
#' Both operands are converted to matrices internally. Supports broadcasting
#' when one operand is a 1x1 scalar matrix.
#'
#' @param e1 A value object or numeric
#' @param e2 A value object or numeric
#'
#' @return A new value object representing the element-wise sum
#' @export
`+.value` <- function(e1, e2)
{
  if (!is_value(e1))
    e1 <- val(e1)
  if (!is_value(e2))
    e2 <- val(e2)

  # Handle broadcasting for different-sized matrices
  bc <- broadcast_matrices(e1$data, e2$data)
  out <- value$new(bc$m1 + bc$m2, list(e1, e2))

  out$backward_fn <- function()
  {
    # Gradient accumulation with broadcasting reduction
    if (is_scalar_matrix(e1$data) && !is_scalar_matrix(out$grad)) {
      e1$grad <- e1$grad + sum(out$grad)
    } else {
      e1$grad <- e1$grad + out$grad
    }
    if (is_scalar_matrix(e2$data) && !is_scalar_matrix(out$grad)) {
      e2$grad <- e2$grad + sum(out$grad)
    } else {
      e2$grad <- e2$grad + out$grad
    }
  }
  out
}

#' Multiplication for value objects
#'
#' Element-wise multiplication of two value objects or a value and numeric.
#' Supports broadcasting when one operand is a 1x1 scalar matrix.
#'
#' @param e1 A value object or numeric
#' @param e2 A value object or numeric
#'
#' @return A new value object representing the element-wise product
#' @export
`*.value` <- function(e1, e2)
{
  if (!is_value(e1))
    e1 <- val(e1)
  if (!is_value(e2))
    e2 <- val(e2)

  # Handle broadcasting
  bc <- broadcast_matrices(e1$data, e2$data)
  out <- value$new(bc$m1 * bc$m2, list(e1, e2))

  out$backward_fn <- function()
  {
    # d(e1*e2)/de1 = e2, d(e1*e2)/de2 = e1
    # With broadcasting, need to reduce gradients for scalar operands
    bc_grad <- broadcast_matrices(e2$data, out$grad)
    g1 <- bc_grad$m1 * bc_grad$m2
    if (is_scalar_matrix(e1$data) && !is_scalar_matrix(g1)) {
      e1$grad <- e1$grad + sum(g1)
    } else {
      e1$grad <- e1$grad + g1
    }

    bc_grad <- broadcast_matrices(e1$data, out$grad)
    g2 <- bc_grad$m1 * bc_grad$m2
    if (is_scalar_matrix(e2$data) && !is_scalar_matrix(g2)) {
      e2$grad <- e2$grad + sum(g2)
    } else {
      e2$grad <- e2$grad + g2
    }
  }
  out
}

#' Power operation for value objects.
#'
#' Element-wise power operation. Supports broadcasting when one operand
#' is a 1x1 scalar matrix.
#'
#' @param b A value object or numeric (base)
#' @param e A value object or numeric (exponent)
#'
#' @return A new value object representing the element-wise power b^e
#' @export
`^.value` <- function(b, e)
{
  if (!is_value(b))
    b <- val(b)

  if (is_value(e))
  {
    bc <- broadcast_matrices(b$data, e$data)
    out <- value$new(bc$m1^bc$m2, list(b, e))
    out$backward_fn <- function()
    {
      # d(b^e)/db = e * b^(e-1)
      bc_be <- broadcast_matrices(b$data, e$data)
      db_local <- bc_be$m2 * bc_be$m1^(bc_be$m2 - 1)
      bc_grad <- broadcast_matrices(db_local, out$grad)
      g_b <- bc_grad$m1 * bc_grad$m2
      if (is_scalar_matrix(b$data) && !is_scalar_matrix(g_b)) {
        b$grad <- b$grad + sum(g_b)
      } else {
        b$grad <- b$grad + g_b
      }

      # d(b^e)/de = b^e * log(b)
      de_local <- bc_be$m1^bc_be$m2 * log(bc_be$m1)
      bc_grad <- broadcast_matrices(de_local, out$grad)
      g_e <- bc_grad$m1 * bc_grad$m2
      if (is_scalar_matrix(e$data) && !is_scalar_matrix(g_e)) {
        e$grad <- e$grad + sum(g_e)
      } else {
        e$grad <- e$grad + g_e
      }
    }
    return(out)
  }
  else
  {
    # e is a plain numeric (scalar or matrix)
    e_mat <- as_matrix(e)
    bc <- broadcast_matrices(b$data, e_mat)
    out <- value$new(bc$m1^bc$m2, list(b))
    out$backward_fn <- function()
    {
      bc_be <- broadcast_matrices(b$data, e_mat)
      db_local <- bc_be$m2 * bc_be$m1^(bc_be$m2 - 1)
      bc_grad <- broadcast_matrices(db_local, out$grad)
      g_b <- bc_grad$m1 * bc_grad$m2
      if (is_scalar_matrix(b$data) && !is_scalar_matrix(g_b)) {
        b$grad <- b$grad + sum(g_b)
      } else {
        b$grad <- b$grad + g_b
      }
    }
    return(out)
  }
}

#' Summation for value objects
#'
#' Sums all elements of value objects, returning a 1x1 value.
#'
#' @param ... value objects and/or numeric values to sum
#' @param na.rm Logical, whether to remove NA values
#'
#' @return A new value object (1x1 matrix) representing the sum
#' @export
sum.value <- function(..., na.rm = FALSE)
{
  values <- list(...)
  if (length(values) == 0)
    return(val(0))

  # Extract value objects and scalars
  value_objs <- Filter(is_value, values)
  scalars <- Filter(Negate(is_value), values)

  # Compute total (sum all elements of all value matrices)
  total <- sum(sapply(value_objs, function(v) sum(v$data, na.rm = na.rm)))
  if (length(scalars) > 0)
    total <- total + sum(unlist(scalars), na.rm = na.rm)

  out <- value$new(total, value_objs)
  out$backward_fn <- function()
  {
    # Gradient flows equally to all elements
    for (child in out$prev)
      child$grad <- child$grad + out$grad[1,1]
  }
  out
}

#' Mean for value objects
#'
#' Computes the arithmetic mean of all elements across value objects.
#'
#' @param x A value object or list of value objects
#' @param ... Additional arguments (for compatibility with base::mean)
#'
#' @return A new value object (1x1 matrix) representing the mean
#' @export
mean.value <- function(x, ...)
{
  # Handle single value object
  if (is_value(x)) {
    n <- length(x$data)
    avg <- sum(x$data) / n
    out <- value$new(avg, list(x))
    out$backward_fn <- function()
    {
      x$grad <- x$grad + out$grad[1,1] / n
    }
    return(out)
  }

  # Handle list of values
  values <- if (is.list(x)) x else list(x)
  value_objs <- Filter(is_value, values)

  if (length(value_objs) == 0)
    return(val(0))

  # Total number of elements
  n <- sum(sapply(value_objs, function(v) length(v$data)))
  total <- sum(sapply(value_objs, function(v) sum(v$data)))
  avg <- total / n

  out <- value$new(avg, value_objs)
  out$backward_fn <- function()
  {
    for (child in out$prev)
      child$grad <- child$grad + out$grad[1,1] / n
  }
  out
}


#' Natural logarithm for value objects
#'
#' Element-wise natural logarithm.
#'
#' @param x A value object
#'
#' @return A new value object representing log(x)
#' @export
log.value <- function(x, ...)
{
  out <- value$new(log(x$data), list(x))
  out$backward_fn <- function()
  {
    x$grad <- x$grad + (1 / x$data) * out$grad
  }
  out
}

#' Exponential function for value objects
#'
#' Element-wise exponential.
#'
#' @param x A value object
#'
#' @return A new value object representing exp(x)
#' @export
exp.value <- function(x)
{
  out <- value$new(exp(x$data), list(x))
  out$backward_fn <- function()
  {
    x$grad <- x$grad + exp(x$data) * out$grad
  }
  out
}

#' Division for value objects
#'
#' Element-wise division. Supports broadcasting when one operand
#' is a 1x1 scalar matrix.
#'
#' @param x A value object or numeric (numerator)
#' @param y A value object or numeric (denominator)
#'
#' @return A new value object representing x / y
#' @export
`/.value` <- function(x, y)
{
  if (!is_value(x))
    x <- val(x)
  if (!is_value(y))
    y <- val(y)

  bc <- broadcast_matrices(x$data, y$data)
  out <- value$new(bc$m1 / bc$m2, list(x, y))

  out$backward_fn <- function()
  {
    bc_xy <- broadcast_matrices(x$data, y$data)

    # d(x/y)/dx = 1/y
    dx_local <- 1 / bc_xy$m2
    bc_grad <- broadcast_matrices(dx_local, out$grad)
    g_x <- bc_grad$m1 * bc_grad$m2
    if (is_scalar_matrix(x$data) && !is_scalar_matrix(g_x)) {
      x$grad <- x$grad + sum(g_x)
    } else {
      x$grad <- x$grad + g_x
    }

    # d(x/y)/dy = -x/y^2
    dy_local <- -bc_xy$m1 / (bc_xy$m2^2)
    bc_grad <- broadcast_matrices(dy_local, out$grad)
    g_y <- bc_grad$m1 * bc_grad$m2
    if (is_scalar_matrix(y$data) && !is_scalar_matrix(g_y)) {
      y$grad <- y$grad + sum(g_y)
    } else {
      y$grad <- y$grad + g_y
    }
  }
  out
}

#' Square root for value objects
#'
#' Element-wise square root.
#'
#' @param x A value object
#'
#' @return A new value object representing sqrt(x)
#' @export
sqrt.value <- function(x)
{
  out <- value$new(sqrt(x$data), list(x))
  out$backward_fn <- function()
  {
    x$grad <- x$grad + (1 / (2 * sqrt(x$data))) * out$grad
  }
  return(out)
}


#' Subtraction for value objects
#'
#' Element-wise subtraction or unary negation. Supports broadcasting
#' when one operand is a 1x1 scalar matrix.
#'
#' @param x A value object or numeric (minuend)
#' @param y A value object or numeric (subtrahend), or NULL for unary negation
#'
#' @return A new value object representing x - y or -x
#' @export
`-.value` <- function(x, y = NULL)
{
  # Handle unary negation
  if (missing(y)) {
    if (!is_value(x))
      x <- val(x)
    out <- value$new(-x$data, list(x))
    out$backward_fn <- function() {
      x$grad <- x$grad - out$grad
    }
    return(out)
  }

  # Handle binary subtraction
  if (!is_value(x))
    x <- val(x)
  if (!is_value(y))
    y <- val(y)

  bc <- broadcast_matrices(x$data, y$data)
  out <- value$new(bc$m1 - bc$m2, list(x, y))

  out$backward_fn <- function()
  {
    if (is_scalar_matrix(x$data) && !is_scalar_matrix(out$grad)) {
      x$grad <- x$grad + sum(out$grad)
    } else {
      x$grad <- x$grad + out$grad
    }
    if (is_scalar_matrix(y$data) && !is_scalar_matrix(out$grad)) {
      y$grad <- y$grad - sum(out$grad)
    } else {
      y$grad <- y$grad - out$grad
    }
  }
  out
}

#' Hyperbolic tangent for value objects
#'
#' Element-wise tanh.
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
    x$grad <- x$grad + (1 - t^2) * out$grad
  }
  out
}

#' Sigmoid activation function for value objects
#'
#' Element-wise sigmoid: 1/(1+exp(-x))
#'
#' @param x A value object
#'
#' @return A new value object representing sigmoid(x)
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
    x$grad <- x$grad + s * (1 - s) * out$grad
  }
  out
}

#' @export
sigmoid.default <- function(x)
{
  1 / (1 + exp(-x))
}

#' ReLU activation function for value objects
#'
#' Element-wise ReLU: max(0, x)
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
  out <- value$new(pmax(0, x$data), list(x))
  out$backward_fn <- function()
  {
    x$grad <- x$grad + (x$data > 0) * out$grad
  }
  out
}

#' @export
relu.default <- function(x)
{
  pmax(0, x)
}

#' Absolute value for value objects
#'
#' Element-wise absolute value.
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
    x$grad <- x$grad + sign(x$data) * out$grad
  }
  out
}

#' Log-gamma function for value objects
#'
#' Element-wise log(gamma(x)).
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
    x$grad <- x$grad + digamma(x$data) * out$grad
  }
  out
}

#' Digamma (psi) function for value objects
#'
#' Element-wise digamma: d/dx log(gamma(x)).
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
    x$grad <- x$grad + trigamma(x$data) * out$grad
  }
  out
}

#' Trigamma function for value objects
#'
#' Element-wise trigamma: second derivative of lgamma.
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
    x$grad <- x$grad + psigamma(x$data, 2) * out$grad
  }
  out
}

#' Log(1+x) for value objects
#'
#' Element-wise log(1+x), numerically stable for small x.
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
    x$grad <- x$grad + (1 / (1 + x$data)) * out$grad
  }
  out
}

#' Sine function for value objects
#'
#' Element-wise sine.
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
    x$grad <- x$grad + cos(x$data) * out$grad
  }
  out
}

#' Cosine function for value objects
#'
#' Element-wise cosine.
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
    x$grad <- x$grad - sin(x$data) * out$grad
  }
  out
}

#' Logit function for value objects
#'
#' Element-wise logit: log(p/(1-p)).
#'
#' @param x A value object representing probabilities
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
    x$grad <- x$grad + (1 / (x$data * (1 - x$data))) * out$grad
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
#' Element-wise softplus: log(1 + exp(x)).
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
  out <- value$new(log1p(exp(x$data)), list(x))
  out$backward_fn <- function()
  {
    # d/dx softplus(x) = sigmoid(x)
    x$grad <- x$grad + (1 / (1 + exp(-x$data))) * out$grad
  }
  out
}

#' @export
softplus.default <- function(x)
{
  log1p(exp(x))
}
