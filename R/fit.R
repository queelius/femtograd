#' Statistical model fitting with automatic differentiation
#'
#' This module provides the `fit()` function for maximum likelihood estimation
#' and the `femtofit` class for representing fitted models. The interface
#' follows standard R conventions, implementing base R generics like
#' `coef()`, `vcov()`, `confint()`, `logLik()`, etc.
#'
#' @name fitting
NULL


#' Fit a model via maximum likelihood
#'
#' Finds the maximum likelihood estimates of parameters and returns a
#' fitted model object with standard errors, confidence intervals, and
#' other inference quantities computed automatically via autodiff.
#'
#' @param loglik A log-likelihood function. Can be specified in two ways:
#'   \itemize{
#'     \item Named arguments: `function(mu, sigma) loglik_normal(mu, sigma, x)`
#'     \item Single parameter argument: `function(p) loglik_normal(p$mu, p$sigma, x)`
#'   }
#'   The function should return a scalar value object.
#' @param params Named numeric vector of initial parameter values,
#'   e.g., `c(mu = 0, sigma = 1)`.
#' @param method Optimization method: "bfgs" (default), "lbfgs", "newton", or "gradient".
#' @param ... Additional arguments passed to the optimizer.
#'
#' @return A `femtofit` object containing:
#'   \itemize{
#'     \item Coefficient estimates (accessible via `coef()`)
#'     \item Variance-covariance matrix (accessible via `vcov()`)
#'     \item Log-likelihood value (accessible via `logLik()`)
#'     \item Standard errors, Hessian, convergence info
#'   }
#'
#' @examples
#' \dontrun{
#' # Generate data
#' set.seed(42)
#' x <- rnorm(100, mean = 5, sd = 2)
#'
#' # Fit using named arguments (recommended)
#' result <- fit(
#'   function(mu, sigma) loglik_normal(mu, sigma, x),
#'   params = c(mu = 0, sigma = 1)
#' )
#'
#' # Or using single parameter argument
#' result <- fit(
#'   function(p) loglik_normal(p$mu, p$sigma, x),
#'   params = c(mu = 0, sigma = 1)
#' )
#'
#' # Standard R generics work
#' coef(result)      # Parameter estimates
#' vcov(result)      # Variance-covariance matrix
#' confint(result)   # 95% confidence intervals
#' logLik(result)    # Log-likelihood (works with AIC/BIC)
#' AIC(result)       # Akaike information criterion
#' summary(result)   # Full summary with p-values
#' }
#'
#' @export
fit <- function(loglik, params, method = c("bfgs", "lbfgs", "newton", "gradient"), ...) {
  method <- match.arg(method)
  call <- match.call()

  # Ensure params are named
  if (is.null(names(params))) {
    names(params) <- paste0("p", seq_along(params))
  }
  param_names <- names(params)
  n <- length(params)

  # Analyze the loglik function signature
  fn_formals <- formals(loglik)
  fn_args <- names(fn_formals)

  # Determine how to call the function
  # Case 1: Function takes named arguments matching param names
  # Case 2: Function takes a single argument (list/params object)
  if (length(fn_args) >= n && all(param_names %in% fn_args)) {
    # Named arguments style: function(mu, sigma, ...)
    wrapped_fn <- function(val_list) {
      args <- setNames(val_list, param_names)
      do.call(loglik, args)
    }
  } else if (length(fn_args) >= 1) {
    # Single argument style: function(p) where p is a named list
    wrapped_fn <- function(val_list) {
      p <- setNames(val_list, param_names)
      loglik(p)
    }
  } else {
    stop("loglik function must accept either named arguments or a single parameter list")
  }

  # Convert initial params to list of value objects
  initial_vals <- lapply(as.list(params), val)

  # Run optimizer
  opt <- switch(method,
    "bfgs" = bfgs(wrapped_fn, initial_vals, maximize = TRUE, ...),
    "lbfgs" = lbfgs(wrapped_fn, initial_vals, maximize = TRUE, ...),
    "newton" = newton_raphson(wrapped_fn, initial_vals, maximize = TRUE, ...),
    "gradient" = gradient_ascent(wrapped_fn, initial_vals, maximize = TRUE, ...)
  )

  # Extract estimates
  estimate <- sapply(opt$params, data)
  names(estimate) <- param_names

  # Compute Hessian if not already computed (newton provides it)
  if (is.null(opt$hessian)) {
    final_vals <- lapply(as.list(estimate), val)
    H <- hessian(wrapped_fn, final_vals)
  } else {
    H <- opt$hessian
  }
  rownames(H) <- colnames(H) <- param_names

  # Compute variance-covariance matrix from Hessian
  # For log-likelihood, vcov = -H^{-1}

  vcov_mat <- tryCatch({
    solve(-H)
  }, error = function(e) {
    warning("Could not invert Hessian; standard errors unavailable")
    matrix(NA_real_, n, n, dimnames = list(param_names, param_names))
  })
  rownames(vcov_mat) <- colnames(vcov_mat) <- param_names

  # Create femtofit object
  femtofit(
    coefficients = estimate,
    vcov = vcov_mat,
    loglik = opt$value,
    hessian = H,
    gradient = opt$gradient,
    converged = opt$converged,
    iterations = opt$iterations,
    call = call
  )
}


#' Constructor for femtofit objects
#'
#' Creates a fitted model object from MLE results. This is primarily
#' an internal constructor; users typically obtain femtofit objects
#' from the `fit()` function.
#'
#' @param coefficients Named numeric vector of parameter estimates.
#' @param vcov Variance-covariance matrix (named).
#' @param loglik Log-likelihood value at the MLE.
#' @param hessian Hessian matrix at the MLE.
#' @param gradient Gradient vector at the MLE (should be near zero).
#' @param converged Logical indicating convergence.
#' @param iterations Number of iterations performed.
#' @param nobs Number of observations (optional).
#' @param call The original function call (optional).
#'
#' @return A `femtofit` object.
#'
#' @export
femtofit <- function(coefficients, vcov, loglik, hessian = NULL,
                     gradient = NULL, converged = TRUE, iterations = NA_integer_,
                     nobs = NULL, call = NULL) {
  structure(
    list(
      coefficients = coefficients,
      vcov = vcov,
      loglik = loglik,
      hessian = hessian,
      gradient = gradient,
      converged = converged,
      iterations = iterations,
      nobs = nobs,
      call = call
    ),
    class = "femtofit"
  )
}


#' @describeIn femtofit Extract coefficient estimates
#' @param object A `femtofit` object.
#' @param ... Additional arguments (ignored).
#' @export
coef.femtofit <- function(object, ...) {
  object$coefficients
}


#' @describeIn femtofit Extract variance-covariance matrix
#' @export
vcov.femtofit <- function(object, ...) {
  object$vcov
}


#' @describeIn femtofit Compute confidence intervals
#' @param parm Parameter names or indices (default: all parameters).
#' @param level Confidence level (default: 0.95).
#' @importFrom stats qnorm
#' @export
confint.femtofit <- function(object, parm, level = 0.95, ...) {
  cf <- coef(object)
  ses <- sqrt(diag(vcov(object)))

  a <- (1 - level) / 2
  z <- qnorm(1 - a)

  ci <- cbind(cf - z * ses, cf + z * ses)
  colnames(ci) <- format_percent(c(a, 1 - a))

  if (!missing(parm)) {
    ci <- ci[parm, , drop = FALSE]
  }

  ci
}


#' @describeIn femtofit Extract log-likelihood (enables AIC/BIC)
#' @export
logLik.femtofit <- function(object, ...) {
  ll <- object$loglik
  attr(ll, "df") <- length(coef(object))
  if (!is.null(object$nobs)) {
    attr(ll, "nobs") <- object$nobs
  }
  class(ll) <- "logLik"
  ll
}


#' @describeIn femtofit Number of observations
#' @importFrom stats nobs
#' @export
nobs.femtofit <- function(object, ...) {
  object$nobs
}


#' @describeIn femtofit Print fitted model
#' @param x A `femtofit` object.
#' @export
print.femtofit <- function(x, ...) {
  cat("femtofit: Maximum Likelihood Estimation\n\n")

  cat("Coefficients:\n")
  print(coef(x))

  cat("\nLog-likelihood:", format(x$loglik, digits = 4), "\n")

  if (!x$converged) {
    cat("Warning: optimization did not converge\n")
  }

  invisible(x)
}


#' @describeIn femtofit Summary of fitted model
#' @export
summary.femtofit <- function(object, ...) {
  cf <- coef(object)
  se <- sqrt(diag(vcov(object)))
  z <- cf / se
  p <- 2 * pnorm(-abs(z))

  coef_table <- cbind(
    Estimate = cf,
    `Std. Error` = se,
    `z value` = z,
    `Pr(>|z|)` = p
  )

  structure(
    list(
      call = object$call,
      coefficients = coef_table,
      loglik = object$loglik,
      aic = AIC(object),
      bic = if (!is.null(object$nobs)) BIC(object) else NA_real_,
      converged = object$converged,
      iterations = object$iterations
    ),
    class = "summary.femtofit"
  )
}


#' @describeIn femtofit Print summary
#' @importFrom stats printCoefmat pnorm AIC BIC
#' @export
print.summary.femtofit <- function(x, ...) {
  cat("femtofit: Maximum Likelihood Estimation\n")

  if (!is.null(x$call)) {
    cat("\nCall:\n")
    print(x$call)
  }

  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE,
               signif.stars = TRUE, na.print = "NA")

  cat("\n---\n")
  cat("Log-likelihood:", format(x$loglik, digits = 4))
  cat("    AIC:", format(x$aic, digits = 4))
  if (!is.na(x$bic)) {
    cat("    BIC:", format(x$bic, digits = 4))
  }
  cat("\n")

  if (x$converged) {
    cat("Converged in", x$iterations, "iterations\n")
  } else {
    cat("Warning: Did not converge after", x$iterations, "iterations\n")
  }

  invisible(x)
}


#' Standard errors from a fitted model
#'
#' Extracts standard errors of coefficient estimates.
#'
#' @param object A `femtofit` object.
#' @param ... Additional arguments (ignored).
#'
#' @return Named numeric vector of standard errors.
#'
#' @export
se <- function(object, ...) {
  UseMethod("se")
}


#' @describeIn se Standard errors for femtofit
#' @export
se.femtofit <- function(object, ...) {
  ses <- sqrt(diag(vcov(object)))
  names(ses) <- names(coef(object))
  ses
}


#' Observed Fisher information matrix
#'
#' Returns the observed Fisher information matrix, which is the
#' negative Hessian of the log-likelihood at the MLE.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments.
#'
#' @return The observed information matrix.
#'
#' @export
observed_info <- function(object, ...) {
  UseMethod("observed_info")
}


#' @describeIn observed_info Observed information for femtofit
#' @export
observed_info.femtofit <- function(object, ...) {
  -object$hessian
}


# Helper function to format percentages for confint column names
format_percent <- function(x) {
  paste0(format(100 * x, trim = TRUE, scientific = FALSE), "%")
}
