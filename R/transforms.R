#' Parameter Transformation Helpers
#'
#' These functions provide clean, differentiable transformations for
#' constrained parameters. They allow optimization in unconstrained space
#' while ensuring parameters stay in valid ranges.
#'
#' @name transforms
NULL


#' Transform to positive values
#'
#' Maps unconstrained values to positive values via exp(x).
#' Use this for parameters like standard deviations, rates, or variances.
#'
#' @param x Unconstrained value (scalar, vector, or value object)
#' @return Positive value (always > 0)
#'
#' @details
#' The transformation is: positive(x) = exp(x)
#'
#' For optimization, work with the unconstrained parameter and transform:
#' \preformatted{
#'   result <- fit(
#'     function(mu, log_sigma) {
#'       sigma <- positive(log_sigma)  # exp(log_sigma)
#'       loglik_normal(mu, sigma, data)
#'     },
#'     params = c(mu = 0, log_sigma = 0)  # log(1) = 0
#'   )
#'   # To recover sigma: exp(coef(result)["log_sigma"])
#' }
#'
#' @examples
#' \dontrun{
#' # Parameter that must be positive
#' log_sigma <- val(-1)
#' sigma <- positive(log_sigma)  # exp(-1) ≈ 0.368
#' data(sigma)
#'
#' # Works in optimization
#' fit(
#'   function(mu, log_sigma) loglik_normal(mu, positive(log_sigma), x),
#'   params = c(mu = 0, log_sigma = 0)
#' )
#' }
#'
#' @export
positive <- function(x) {
  exp(x)
}


#' Transform to probability values
#'
#' Maps unconstrained values to (0, 1) via the sigmoid function.
#' Use this for probability parameters.
#'
#' @param x Unconstrained value (scalar, vector, or value object)
#' @return Value in (0, 1)
#'
#' @details
#' The transformation is: probability(x) = 1 / (1 + exp(-x)) = sigmoid(x)
#'
#' For optimization:
#' \preformatted{
#'   result <- fit(
#'     function(logit_p) {
#'       p <- probability(logit_p)
#'       loglik_bernoulli(p, data)
#'     },
#'     params = c(logit_p = 0)  # sigmoid(0) = 0.5
#'   )
#'   # To recover p: sigmoid(coef(result)["logit_p"])
#' }
#'
#' @examples
#' \dontrun{
#' # Parameter that must be in (0, 1)
#' logit_p <- val(2)
#' p <- probability(logit_p)  # sigmoid(2) ≈ 0.88
#' data(p)
#' }
#'
#' @export
probability <- function(x) {
  sigmoid(x)
}


#' Transform to bounded interval
#'
#' Maps unconstrained values to a bounded interval (lower, upper).
#' Use this for parameters with both lower and upper bounds.
#'
#' @param x Unconstrained value (scalar, vector, or value object)
#' @param lower Lower bound of the interval
#' @param upper Upper bound of the interval
#' @return Value in (lower, upper)
#'
#' @details
#' The transformation is:
#' \preformatted{
#'   bounded(x, a, b) = a + (b - a) * sigmoid(x)
#' }
#'
#' This maps (-∞, ∞) to (a, b) smoothly.
#'
#' @examples
#' \dontrun{
#' # Parameter in (0, 10)
#' raw <- val(0)
#' param <- bounded(raw, 0, 10)  # sigmoid(0) * 10 = 5
#' data(param)
#'
#' # Correlation parameter in (-1, 1)
#' raw_rho <- val(1)
#' rho <- bounded(raw_rho, -1, 1)
#' }
#'
#' @export
bounded <- function(x, lower, upper) {
  stopifnot(lower < upper)
  range <- upper - lower
  lower + range * sigmoid(x)
}


#' Transform to lower-bounded interval
#'
#' Maps unconstrained values to (lower, ∞).
#' Equivalent to lower + positive(x).
#'
#' @param x Unconstrained value
#' @param lower Lower bound
#' @return Value > lower
#'
#' @details
#' The transformation is: lower_bounded(x, a) = a + exp(x)
#'
#' @examples
#' \dontrun{
#' # Parameter > 2
#' raw <- val(0)
#' param <- lower_bounded(raw, 2)  # 2 + exp(0) = 3
#' }
#'
#' @export
lower_bounded <- function(x, lower) {

  lower + exp(x)
}


#' Transform to upper-bounded interval
#'
#' Maps unconstrained values to (-∞, upper).
#' Equivalent to upper - positive(x).
#'
#' @param x Unconstrained value
#' @param upper Upper bound
#' @return Value < upper
#'
#' @details
#' The transformation is: upper_bounded(x, b) = b - exp(x)
#'
#' @examples
#' \dontrun{
#' # Parameter < 10
#' raw <- val(0)
#' param <- upper_bounded(raw, 10)  # 10 - exp(0) = 9
#' }
#'
#' @export
upper_bounded <- function(x, upper) {
  upper - exp(x)
}


#' Inverse transforms for recovering original scale
#'
#' These functions convert fitted unconstrained parameters back to
#' their natural scale.
#'
#' @name inverse_transforms
NULL


#' Inverse of positive transform
#'
#' Converts a positive value back to unconstrained scale.
#'
#' @param x Positive value
#' @return Unconstrained value (log of x)
#'
#' @export
inv_positive <- function(x) {
  log(x)
}


#' Inverse of probability transform
#'
#' Converts a probability to unconstrained scale (logit).
#'
#' @param x Value in (0, 1)
#' @return Unconstrained value (logit of x)
#'
#' @export
inv_probability <- function(x) {
  logit(x)
}


#' Inverse of bounded transform
#'
#' Converts a bounded value back to unconstrained scale.
#'
#' @param x Value in (lower, upper)
#' @param lower Lower bound
#' @param upper Upper bound
#' @return Unconstrained value
#'
#' @export
inv_bounded <- function(x, lower, upper) {
  stopifnot(lower < upper)
  # x = lower + (upper - lower) * sigmoid(raw)

# sigmoid(raw) = (x - lower) / (upper - lower)
  # raw = logit((x - lower) / (upper - lower))
  p <- (x - lower) / (upper - lower)
  logit(p)
}
