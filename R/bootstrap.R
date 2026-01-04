#' Bootstrap Inference
#'
#' Functions for bootstrap-based inference on fitted models.
#' Bootstrap methods provide non-parametric standard errors and
#' confidence intervals without relying on asymptotic normality.
#'
#' @name bootstrap
NULL


#' Bootstrap standard errors and confidence intervals
#'
#' Performs bootstrap resampling to estimate the sampling distribution
#' of parameter estimates.
#'
#' @param loglik_maker A function that takes data and returns a log-likelihood
#'   function suitable for `fit()`. See examples.
#' @param data The data (vector, matrix, or data frame).
#' @param params Initial parameter values for `fit()`.
#' @param n_boot Number of bootstrap replications (default 500).
#' @param progress Logical; show progress messages (default TRUE).
#' @param ... Additional arguments passed to `fit()`.
#'
#' @return A `bootstrap_result` object containing:
#'   \describe{
#'     \item{estimates}{Matrix of bootstrap estimates (n_boot x n_params)}
#'     \item{original}{Original parameter estimates}
#'     \item{se}{Bootstrap standard errors}
#'     \item{bias}{Estimated bias}
#'   }
#'
#' @details
#' The `loglik_maker` function should accept data and return a log-likelihood
#' function that can be passed to `fit()`. This design allows bootstrap
#' to work with any model specification.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' x <- rnorm(50, mean = 5, sd = 2)
#'
#' # Define a function that creates the loglik function from data
#' make_loglik <- function(data) {
#'   function(mu, log_sigma) {
#'     loglik_normal(mu, exp(log_sigma), data)
#'   }
#' }
#'
#' # Bootstrap
#' boot <- bootstrap_fit(
#'   make_loglik,
#'   data = x,
#'   params = c(mu = 0, log_sigma = 0),
#'   n_boot = 200
#' )
#'
#' # Results
#' boot$se          # Bootstrap SEs
#' confint(boot)    # Bootstrap CIs
#'
#' # Compare with Wald SEs
#' result <- fit(make_loglik(x), params = c(mu = 0, log_sigma = 0))
#' se(result)       # Wald SEs
#' }
#'
#' @export
bootstrap_fit <- function(loglik_maker, data, params, n_boot = 500,
                          progress = TRUE, ...) {

  # Fit original model
  loglik_fn <- loglik_maker(data)
  original_fit <- fit(loglik_fn, params = params, ...)

  if (!original_fit$converged) {
    warning("Original model did not converge")
  }

  original <- coef(original_fit)
  n_params <- length(original)
  param_names <- names(original)

  # Determine sample size
  if (is.vector(data)) {
    n <- length(data)
    resample_fn <- function() sample(data, n, replace = TRUE)
  } else if (is.matrix(data) || is.data.frame(data)) {
    n <- nrow(data)
    resample_fn <- function() {
      idx <- sample(n, n, replace = TRUE)
      data[idx, , drop = FALSE]
    }
  } else {
    stop("data must be a vector, matrix, or data frame")
  }

  # Storage for bootstrap estimates
  boot_estimates <- matrix(NA_real_, nrow = n_boot, ncol = n_params)
  colnames(boot_estimates) <- param_names

  n_failed <- 0

  # Bootstrap loop
  for (b in seq_len(n_boot)) {
    if (progress && b %% 100 == 0) {
      message(sprintf("Bootstrap: %d/%d", b, n_boot))
    }

    boot_data <- resample_fn()

    boot_fit <- tryCatch({
      boot_loglik <- loglik_maker(boot_data)
      fit(boot_loglik, params = original, ...)
    }, error = function(e) {
      NULL
    })

    if (!is.null(boot_fit) && boot_fit$converged) {
      boot_estimates[b, ] <- coef(boot_fit)
    } else {
      n_failed <- n_failed + 1
    }
  }

  # Remove failed fits
  valid_rows <- !apply(is.na(boot_estimates), 1, any)
  boot_estimates <- boot_estimates[valid_rows, , drop = FALSE]
  n_successful <- sum(valid_rows)

  if (n_successful < 50) {
    warning(sprintf("Only %d successful bootstrap samples. Results may be unreliable.",
                    n_successful))
  }

  # Compute bootstrap statistics
  boot_means <- colMeans(boot_estimates)
  boot_se <- apply(boot_estimates, 2, sd)
  boot_bias <- boot_means - original

  structure(
    list(
      estimates = boot_estimates,
      original = original,
      se = boot_se,
      bias = boot_bias,
      mean = boot_means,
      n_boot = n_boot,
      n_successful = n_successful,
      n_failed = n_failed
    ),
    class = "bootstrap_result"
  )
}


#' Confidence intervals from bootstrap
#'
#' Computes confidence intervals using bootstrap methods.
#'
#' @param object A `bootstrap_result` object.
#' @param parm Parameter names or indices (default: all).
#' @param level Confidence level (default 0.95).
#' @param type Type of interval: "percentile" (default), "basic", or "normal".
#' @param ... Additional arguments (ignored).
#'
#' @return Matrix of confidence intervals.
#'
#' @details
#' Available methods:
#' \describe{
#'   \item{percentile}{Uses quantiles of bootstrap distribution directly}
#'   \item{basic}{Uses 2*theta_hat - quantiles (pivot method)}
#'   \item{normal}{Uses bootstrap SE with normal quantiles}
#' }
#'
#' @importFrom stats quantile qnorm
#' @export
confint.bootstrap_result <- function(object, parm = NULL, level = 0.95,
                                     type = c("percentile", "basic", "normal"),
                                     ...) {
  type <- match.arg(type)

  if (is.null(parm)) {
    parm <- colnames(object$estimates)
  } else if (is.numeric(parm)) {
    parm <- colnames(object$estimates)[parm]
  }

  alpha <- (1 - level) / 2
  probs <- c(alpha, 1 - alpha)

  ci <- matrix(NA_real_, nrow = length(parm), ncol = 2)
  rownames(ci) <- parm
  colnames(ci) <- paste0(format(100 * probs, trim = TRUE), "%")

  for (i in seq_along(parm)) {
    p <- parm[i]
    boot_vals <- object$estimates[, p]
    theta_hat <- object$original[p]

    if (type == "percentile") {
      ci[i, ] <- quantile(boot_vals, probs = probs, na.rm = TRUE)
    } else if (type == "basic") {
      q <- quantile(boot_vals, probs = rev(probs), na.rm = TRUE)
      ci[i, ] <- 2 * theta_hat - q
    } else if (type == "normal") {
      z <- qnorm(1 - alpha)
      se <- object$se[p]
      ci[i, ] <- theta_hat + c(-z, z) * se
    }
  }

  ci
}


#' Print method for bootstrap results
#'
#' @param x A `bootstrap_result` object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.bootstrap_result <- function(x, ...) {
  cat("Bootstrap Results\n")
  cat("-----------------\n")
  cat("Bootstrap samples:", x$n_successful, "of", x$n_boot,
      sprintf("(%.1f%% success)\n", 100 * x$n_successful / x$n_boot))
  cat("\n")

  # Create summary table
  tab <- data.frame(
    Original = x$original,
    Boot.SE = x$se,
    Bias = x$bias
  )

  print(tab, digits = 4)

  invisible(x)
}


#' Summary for bootstrap results
#'
#' @param object A `bootstrap_result` object.
#' @param level Confidence level (default 0.95).
#' @param ... Additional arguments (ignored).
#'
#' @export
summary.bootstrap_result <- function(object, level = 0.95, ...) {
  ci <- confint(object, level = level)

  tab <- data.frame(
    Estimate = object$original,
    Boot.SE = object$se,
    Lower = ci[, 1],
    Upper = ci[, 2]
  )

  structure(
    list(
      coefficients = tab,
      level = level,
      n_boot = object$n_boot,
      n_successful = object$n_successful
    ),
    class = "summary.bootstrap_result"
  )
}


#' @export
print.summary.bootstrap_result <- function(x, ...) {
  cat("Bootstrap Summary\n")
  cat("-----------------\n")
  cat("Successful samples:", x$n_successful, "of", x$n_boot, "\n")
  cat(sprintf("Confidence level: %.0f%%\n\n", 100 * x$level))

  print(x$coefficients, digits = 4)

  invisible(x)
}
