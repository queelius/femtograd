#' Numerical stability utilities for automatic differentiation
#'
#' These functions provide numerically stable implementations of common
#' operations that can overflow or underflow with naive implementations.
#'
#' @name stability
NULL


#' Log-Sum-Exp (numerically stable)
#'
#' Computes log(sum(exp(x))) in a numerically stable way by factoring out
#' the maximum value: log(sum(exp(x))) = m + log(sum(exp(x - m)))
#' where m = max(x).
#'
#' @param ... value objects and/or numeric vectors to include in the sum
#' @param na.rm Logical, whether to remove NA values
#'
#' @return A value object representing log(sum(exp(x)))
#'
#' @details
#' The naive computation of log(sum(exp(x))) overflows when x contains
#' large values (exp(710) = Inf in double precision) and underflows when
#' x contains very negative values. This implementation is stable for
#' any finite input.
#'
#' @examples
#' \dontrun{
#' x <- val(c(1000, 1001, 1002))  # Would overflow with naive exp()
#' logsumexp(x)  # Returns ~1002.41 correctly
#' }
#'
#' @export
logsumexp <- function(..., na.rm = FALSE)
{
  UseMethod("logsumexp")
}

#' @export
logsumexp.value <- function(..., na.rm = FALSE)
{
  args <- list(...)

  # Collect all data values to find the max
  all_data <- unlist(lapply(args, function(a) {
    if (is_value(a)) a$data else a
  }))

  if (length(all_data) == 0) return(val(-Inf))

  m <- max(all_data, na.rm = na.rm)

  # Handle edge case: all -Inf
  if (is.infinite(m) && m < 0) return(val(-Inf))

  # Compute exp(x - m) for each argument, then sum all elements
  # Use sum() to reduce vectors to scalars
  exp_terms <- lapply(args, function(a) {
    if (is_value(a)) {
      exp(a - m)
    } else {
      val(sum(exp(a - m), na.rm = na.rm))
    }
  })

  # Sum all terms - use sum() on each to handle vector-valued value objects
  # Then combine with addition
  if (length(exp_terms) == 1) {
    total <- sum(exp_terms[[1]])
  } else {
    # Sum each term individually, then add them together
    summed <- lapply(exp_terms, function(e) {
      if (is_value(e) && length(e$data) > 1) sum(e) else e
    })
    total <- Reduce(`+`, summed)
  }

  # Return m + log(total)
  m + log(total)
}

#' @export
logsumexp.numeric <- function(..., na.rm = FALSE)
{
  args <- list(...)
  all_vals <- unlist(args)

  if (length(all_vals) == 0) return(-Inf)

  m <- max(all_vals, na.rm = na.rm)
  if (is.infinite(m) && m < 0) return(-Inf)

  m + log(sum(exp(all_vals - m), na.rm = na.rm))
}

#' @export
logsumexp.default <- function(..., na.rm = FALSE)
{
  logsumexp.numeric(..., na.rm = na.rm)
}


#' Softmax function (numerically stable)
#'
#' Computes softmax(x)_i = exp(x_i) / sum(exp(x)) in a numerically stable way.
#'
#' @param x A value object or numeric vector
#'
#' @return A value object (or numeric) with softmax probabilities
#'
#' @details
#' Uses the identity softmax(x) = softmax(x - max(x)) to prevent overflow.
#'
#' @export
softmax <- function(x)
{
  UseMethod("softmax")
}

#' @export
softmax.value <- function(x)
{
  m <- max(x$data)
  shifted <- x - m
  e <- exp(shifted)
  e / sum(e)
}

#' @export
softmax.numeric <- function(x)
{
  m <- max(x)
  e <- exp(x - m)
  e / sum(e)
}

#' @export
softmax.default <- function(x)
{
  softmax.numeric(x)
}


#' Safe logarithm (handles zeros)
#'
#' Computes log(x) with protection against log(0) = -Inf.
#' Optionally clamps input to a minimum value.
#'
#' @param x A value object or numeric
#' @param eps Minimum value to clamp x to (default: .Machine$double.eps)
#'
#' @return log(max(x, eps))
#'
#' @details
#' This is useful when computing log-likelihoods where probabilities
#' might numerically become zero. The gradient is computed as if
#' the clamping didn't happen (straight-through estimator) when x > eps.
#'
#' @export
log_safe <- function(x, eps = .Machine$double.eps)
{
  UseMethod("log_safe")
}

#' @export
log_safe.value <- function(x, eps = .Machine$double.eps)
{
  # Clamp data but compute gradient as normal where x > eps
  clamped_data <- pmax(x$data, eps)
  out <- value$new(log(clamped_data), list(x))
  out$backward_fn <- function()
  {
    # Gradient is 1/x where x > eps, 1/eps where x <= eps
    safe_x <- pmax(x$data, eps)
    x$grad <<- x$grad + (1 / safe_x) * out$grad
  }
  out
}

#' @export
log_safe.numeric <- function(x, eps = .Machine$double.eps)
{
  log(pmax(x, eps))
}

#' @export
log_safe.default <- function(x, eps = .Machine$double.eps)
{
  log_safe.numeric(x, eps)
}


#' Safe division (handles division by zero)
#'
#' Computes x / y with protection against division by zero.
#'
#' @param x Numerator (value object or numeric)
#' @param y Denominator (value object or numeric)
#' @param eps Minimum absolute value for denominator (default: .Machine$double.eps)
#'
#' @return x / sign(y) * max(abs(y), eps)
#'
#' @export
div_safe <- function(x, y, eps = .Machine$double.eps)
{
  UseMethod("div_safe")
}

#' @export
div_safe.value <- function(x, y, eps = .Machine$double.eps)
{
  if (!is_value(x)) x <- val(x)
  if (!is_value(y)) y <- val(y)

  # Safe denominator: preserve sign, clamp magnitude
  y_safe_data <- sign(y$data) * pmax(abs(y$data), eps)
  y_safe_data[y_safe_data == 0] <- eps  # Handle exact zero

  out <- value$new(x$data / y_safe_data, list(x, y))
  out$backward_fn <- function()
  {
    y_safe <- sign(y$data) * pmax(abs(y$data), eps)
    y_safe[y_safe == 0] <- eps
    x$grad <<- x$grad + (1 / y_safe) * out$grad
    y$grad <<- y$grad - (x$data / (y_safe^2)) * out$grad
  }
  out
}

#' @export
div_safe.numeric <- function(x, y, eps = .Machine$double.eps)
{
  y_safe <- sign(y) * pmax(abs(y), eps)
  y_safe[y_safe == 0] <- eps
  x / y_safe
}

#' @export
div_safe.default <- function(x, y, eps = .Machine$double.eps)
{
  div_safe.numeric(x, y, eps)
}


#' Stable sigmoid function
#'
#' Computes sigmoid(x) = 1/(1+exp(-x)) with overflow protection.
#'
#' @param x A value object or numeric
#'
#' @return Sigmoid values in (0, 1)
#'
#' @details
#' For large positive x, exp(-x) underflows to 0, giving sigmoid = 1 (correct).
#' For large negative x, we use sigmoid(x) = exp(x)/(1+exp(x)) to avoid overflow.
#'
#' @export
sigmoid_stable <- function(x)
{
  UseMethod("sigmoid_stable")
}

#' @export
sigmoid_stable.value <- function(x)
{
  # Use stable computation based on sign of x
  # For x >= 0: 1 / (1 + exp(-x))
  # For x < 0: exp(x) / (1 + exp(x))
  pos_mask <- x$data >= 0

  result_data <- ifelse(
    pos_mask,
    1 / (1 + exp(-x$data)),
    exp(x$data) / (1 + exp(x$data))
  )

  out <- value$new(result_data, list(x))
  out$backward_fn <- function()
  {
    s <- out$data
    x$grad <<- x$grad + s * (1 - s) * out$grad
  }
  out
}

#' @export
sigmoid_stable.numeric <- function(x)
{
  ifelse(x >= 0,
         1 / (1 + exp(-x)),
         exp(x) / (1 + exp(x)))
}

#' @export
sigmoid_stable.default <- function(x)
{
  sigmoid_stable.numeric(x)
}


#' Stable exp function (with overflow protection)
#'
#' Computes exp(x) with optional clamping to prevent Inf.
#'
#' @param x A value object or numeric
#' @param max_val Maximum value for x to prevent overflow (default: 709)
#'
#' @return exp(min(x, max_val))
#'
#' @details
#' In double precision, exp(710) overflows to Inf. This function clamps
#' the input to prevent overflow while maintaining correct gradients.
#'
#' @export
exp_safe <- function(x, max_val = 709)
{
  UseMethod("exp_safe")
}

#' @export
exp_safe.value <- function(x, max_val = 709)
{
  clamped_data <- pmin(x$data, max_val)
  result <- exp(clamped_data)

  out <- value$new(result, list(x))
  out$backward_fn <- function()
  {
    # Gradient is exp(x) where x <= max_val, exp(max_val) otherwise
    x$grad <<- x$grad + result * out$grad
  }
  out
}

#' @export
exp_safe.numeric <- function(x, max_val = 709)
{
  exp(pmin(x, max_val))
}

#' @export
exp_safe.default <- function(x, max_val = 709)
{
  exp_safe.numeric(x, max_val)
}


#' Log1p with underflow protection
#'
#' Computes log(1 + x) with protection for very small x.
#'
#' @param x A value object or numeric
#'
#' @return log(1 + x), accurate even for x near zero
#'
#' @details
#' R's built-in log1p is already stable, but this version adds
#' autodiff support and handles the case where 1 + x might underflow
#' to 1 for very small x.
#'
#' @export
log1p_safe <- function(x)
{
  UseMethod("log1p_safe")
}

#' @export
log1p_safe.value <- function(x)
{
  # R's log1p is already numerically stable
  log1p(x)
}

#' @export
log1p_safe.numeric <- function(x)
{
  log1p(x)
}

#' @export
log1p_safe.default <- function(x)
{
  log1p(x)
}


#' Log-sigmoid (numerically stable)
#'
#' Computes log(sigmoid(x)) = -log(1 + exp(-x)) stably.
#'
#' @param x A value object or numeric
#'
#' @return log(sigmoid(x))
#'
#' @details
#' Direct computation of log(sigmoid(x)) fails for large negative x
#' (sigmoid underflows to 0). This uses:
#' - For x >= 0: -log(1 + exp(-x)) = -softplus(-x)
#' - For x < 0: x - log(1 + exp(x)) = x - softplus(x)
#'
#' @export
log_sigmoid <- function(x)
{
  UseMethod("log_sigmoid")
}

#' @export
log_sigmoid.value <- function(x)
{
  # log(sigmoid(x)) = -softplus(-x)
  -softplus(-x)
}

#' @export
log_sigmoid.numeric <- function(x)
{
  # Stable computation
  ifelse(x >= 0,
         -log1p(exp(-x)),
         x - log1p(exp(x)))
}

#' @export
log_sigmoid.default <- function(x)
{
  log_sigmoid.numeric(x)
}
