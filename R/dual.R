#' dual R6 class for forward-mode automatic differentiation
#'
#' Represents a dual number (primal, tangent) for computing directional
#' derivatives via forward-mode AD. When combined with reverse-mode (value),
#' enables Hessian computation via forward-over-reverse.
#'
#' @field primal The function value (can be numeric or value object)
#' @field tangent The directional derivative (can be numeric or value object)
#'
#' @details
#' Forward-mode AD propagates derivatives alongside values during the forward
#' pass. For a function f(x), setting tangent(x) = 1 gives f'(x) in the
#' tangent of the output.
#'
#' For Hessian computation (forward-over-reverse):
#' - primal is a \code{value} object (for reverse-mode gradient)
#' - tangent is also a \code{value} object (tracks d/dx of the gradient)
#' - After backward() on primal, tangent of grad gives Hessian entries
#'
#' @import R6
#' @export
dual <- R6Class("dual",
  public = list(
    primal = NULL,
    tangent = NULL,

    #' @description Create a new dual number
    #' @param primal The primal value (numeric or value object)
    #' @param tangent The tangent (derivative direction), default 0
    initialize = function(primal, tangent = 0)
    {
      self$primal <- primal
      self$tangent <- tangent
    }
  )
)

#' Create a dual number
#'
#' @param primal The primal value
#' @param tangent The tangent (directional derivative), default 0
#'
#' @return A dual object
#' @export
dual_num <- function(primal, tangent = 0)
{
  dual$new(primal, tangent)
}

#' Check if object is a dual number
#' @param x Object to check
#' @return TRUE if x is a dual object
#' @export
is_dual <- function(x) inherits(x, "dual")

#' Extract primal from dual or return value unchanged
#' @param x A dual or non-dual object
#' @return The primal value
#' @export
primal <- function(x)
{
  if (is_dual(x)) x$primal else x
}

#' Extract tangent from dual or return 0
#' @param x A dual or non-dual object
#' @return The tangent value
#' @export
tangent <- function(x)
{
  if (is_dual(x)) x$tangent else 0
}

# Helper to check if tangent is a scalar zero (safe for vectors)
is_zero_tangent <- function(x)
{
  identical(x, 0) || identical(x, 0L)
}

# Arithmetic operations for dual numbers
# Rule: (a, a') op (b, b') where a' and b' are tangents

#' @export
`+.dual` <- function(e1, e2)
{
  if (!is_dual(e2)) e2 <- dual_num(e2, 0)
  if (!is_dual(e1)) e1 <- dual_num(e1, 0)
  a <- primal(e1)
  b <- primal(e2)
  da <- tangent(e1)
  db <- tangent(e2)
  # (a + b, a' + b')
  primal_result <- a + b
  # Handle zero tangents carefully
  if (is_zero_tangent(da)) {
    tangent_result <- db
  } else if (is_zero_tangent(db)) {
    tangent_result <- da
  } else {
    tangent_result <- da + db
  }
  dual_num(primal_result, tangent_result)
}

#' @export
`-.dual` <- function(e1, e2 = NULL)
{
  if (is.null(e2)) {
    # Unary negation
    da <- tangent(e1)
    neg_tangent <- if (is_zero_tangent(da)) 0 else -da
    return(dual_num(-primal(e1), neg_tangent))
  }
  if (!is_dual(e2)) e2 <- dual_num(e2, 0)
  if (!is_dual(e1)) e1 <- dual_num(e1, 0)
  a <- primal(e1)
  b <- primal(e2)
  da <- tangent(e1)
  db <- tangent(e2)
  # (a - b, a' - b')
  primal_result <- a - b
  if (is_zero_tangent(da)) {
    tangent_result <- if (is_zero_tangent(db)) 0 else -db
  } else if (is_zero_tangent(db)) {
    tangent_result <- da
  } else {
    tangent_result <- da - db
  }
  dual_num(primal_result, tangent_result)
}

#' @export
`*.dual` <- function(e1, e2)
{
  if (!is_dual(e2)) e2 <- dual_num(e2, 0)
  if (!is_dual(e1)) e1 <- dual_num(e1, 0)
  a <- primal(e1)
  b <- primal(e2)
  da <- tangent(e1)
  db <- tangent(e2)
  # (a * b, a' * b + a * b')
  # Ensure value objects come first for proper dispatch
  primal_result <- a * b
  term1 <- if (is_zero_tangent(da)) 0 else da * b
  term2 <- if (is_zero_tangent(db)) 0 else a * db
  tangent_result <- term1 + term2
  dual_num(primal_result, tangent_result)
}

#' @export
`/.dual` <- function(e1, e2)
{
  if (!is_dual(e2)) e2 <- dual_num(e2, 0)
  if (!is_dual(e1)) e1 <- dual_num(e1, 0)
  a <- primal(e1)
  b <- primal(e2)
  da <- tangent(e1)
  db <- tangent(e2)
  # (a/b, (a'*b - a*b') / b^2)
  primal_result <- a / b
  b2 <- b * b
  term1 <- if (is_zero_tangent(da)) 0 else da * b
  term2 <- if (is_zero_tangent(db)) 0 else a * db
  tangent_result <- (term1 - term2) / b2
  dual_num(primal_result, tangent_result)
}

#' @export
`^.dual` <- function(e1, e2)
{
  if (!is_dual(e2)) e2 <- dual_num(e2, 0)
  if (!is_dual(e1)) e1 <- dual_num(e1, 0)
  a <- primal(e1)
  b <- primal(e2)
  da <- tangent(e1)
  db <- tangent(e2)
  # d/dx (a^b) = b*a^(b-1)*a' + a^b*log(a)*b'
  # Note: when a,da are value objects and b,db are scalars,
  # we need value to come first for correct dispatch
  primal_result <- a^b
  # Compute tangent: da * (b * a^(b-1)) + db * (a^b * log(a))
  term1 <- if (is_zero_tangent(da)) 0 else da * (a^(b - 1) * b)
  term2 <- if (is_zero_tangent(db)) 0 else (a^b * log(a)) * db
  tangent_result <- term1 + term2
  dual_num(primal_result, tangent_result)
}

#' @export
log.dual <- function(x)
{
  a <- primal(x)
  da <- tangent(x)
  # d/dx log(a) = a'/a
  primal_result <- log(a)
  tangent_result <- if (is_zero_tangent(da)) 0 else da / a
  dual_num(primal_result, tangent_result)
}

#' @export
exp.dual <- function(x)
{
  a <- primal(x)
  da <- tangent(x)
  # d/dx exp(a) = exp(a) * a'
  e <- exp(a)
  tangent_result <- if (is_zero_tangent(da)) 0 else e * da
  dual_num(e, tangent_result)
}

#' @export
sqrt.dual <- function(x)
{
  a <- primal(x)
  da <- tangent(x)
  # d/dx sqrt(a) = a' / (2*sqrt(a))
  s <- sqrt(a)
  tangent_result <- if (is_zero_tangent(da)) 0 else da / (s * 2)
  dual_num(s, tangent_result)
}

#' @export
sin.dual <- function(x)
{
  a <- primal(x)
  da <- tangent(x)
  tangent_result <- if (is_zero_tangent(da)) 0 else cos(a) * da
  dual_num(sin(a), tangent_result)
}

#' @export
cos.dual <- function(x)
{
  a <- primal(x)
  da <- tangent(x)
  tangent_result <- if (is_zero_tangent(da)) 0 else -sin(a) * da
  dual_num(cos(a), tangent_result)
}

#' @export
tanh.dual <- function(x)
{
  a <- primal(x)
  da <- tangent(x)
  t <- tanh(a)
  if (is_value(a)) {
    tangent_result <- if (is_zero_tangent(da)) 0 else (val(1) - t * t) * da
  } else {
    tangent_result <- if (is_zero_tangent(da)) 0 else (1 - t * t) * da
  }
  dual_num(t, tangent_result)
}

#' @export
abs.dual <- function(x)
{
  a <- primal(x)
  da <- tangent(x)
  tangent_result <- if (is_zero_tangent(da)) 0 else sign(a) * da
  dual_num(abs(a), tangent_result)
}

#' @export
lgamma.dual <- function(x)
{
  a <- primal(x)
  da <- tangent(x)
  tangent_result <- if (is_zero_tangent(da)) 0 else digamma(a) * da
  dual_num(lgamma(a), tangent_result)
}

#' @export
digamma.dual <- function(x)
{
  a <- primal(x)
  da <- tangent(x)
  tangent_result <- if (is_zero_tangent(da)) 0 else trigamma(a) * da
  dual_num(digamma(a), tangent_result)
}

#' @export
log1p.dual <- function(x)
{
  a <- primal(x)
  da <- tangent(x)
  if (is_value(a)) {
    tangent_result <- if (is_zero_tangent(da)) 0 else da / (val(1) + a)
  } else {
    tangent_result <- if (is_zero_tangent(da)) 0 else da / (1 + a)
  }
  dual_num(log1p(a), tangent_result)
}

#' @export
sigmoid.dual <- function(x)
{
  a <- primal(x)
  da <- tangent(x)
  # Use value wrapper only if primal is a value object
  if (is_value(a)) {
    s <- val(1) / (val(1) + exp(-a))
    tangent_result <- if (is_zero_tangent(da)) 0 else s * (val(1) - s) * da
  } else {
    s <- 1 / (1 + exp(-a))
    tangent_result <- if (is_zero_tangent(da)) 0 else s * (1 - s) * da
  }
  dual_num(s, tangent_result)
}

#' @export
relu.dual <- function(x)
{
  a <- primal(x)
  da <- tangent(x)
  # For value objects, need to extract numeric for comparison
  a_num <- if (is_value(a)) data(a) else a
  primal_result <- if (a_num > 0) a else val(0)
  tangent_result <- if (a_num > 0) da else 0
  dual_num(primal_result, tangent_result)
}

#' @export
softplus.dual <- function(x)
{
  a <- primal(x)
  da <- tangent(x)
  if (is_value(a)) {
    s <- val(1) / (val(1) + exp(-a))  # sigmoid
    tangent_result <- if (is_zero_tangent(da)) 0 else s * da
  } else {
    s <- 1 / (1 + exp(-a))  # sigmoid
    tangent_result <- if (is_zero_tangent(da)) 0 else s * da
  }
  dual_num(log1p(exp(a)), tangent_result)
}

#' @export
logit.dual <- function(x)
{
  a <- primal(x)
  da <- tangent(x)
  if (is_value(a)) {
    tangent_result <- if (is_zero_tangent(da)) 0 else da / (a * (val(1) - a))
    primal_result <- log(a / (val(1) - a))
  } else {
    tangent_result <- if (is_zero_tangent(da)) 0 else da / (a * (1 - a))
    primal_result <- log(a / (1 - a))
  }
  dual_num(primal_result, tangent_result)
}

#' @export
print.dual <- function(x, ...)
{
  cat("dual(primal =", format(primal(x)), ", tangent =", format(tangent(x)), ")\n")
}

#' Sum for dual numbers
#' @export
sum.dual <- function(..., na.rm = FALSE)
{
  args <- list(...)
  total_primal <- 0
  total_tangent <- 0
  for (arg in args) {
    if (is_dual(arg)) {
      # Use sum() to reduce vector primals/tangents to scalars
      total_primal <- total_primal + sum(primal(arg), na.rm = na.rm)
      total_tangent <- total_tangent + sum(tangent(arg), na.rm = na.rm)
    } else {
      total_primal <- total_primal + sum(arg, na.rm = na.rm)
    }
  }
  dual_num(total_primal, total_tangent)
}
