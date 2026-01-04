#' Model Diagnostics
#'
#' Functions for checking model fit quality, convergence, and
#' potential numerical issues.
#'
#' @name diagnostics
NULL


#' Check Hessian properties
#'
#' Examines the Hessian matrix for numerical issues that may indicate
#' problems with the optimization or model specification.
#'
#' @param object A `femtofit` object.
#'
#' @return A `hessian_check` object containing:
#'   \describe{
#'     \item{is_negative_definite}{Logical: is -H positive definite? (required for MLE)}
#'     \item{eigenvalues}{Eigenvalues of -H (observed information)}
#'     \item{condition_number}{Ratio of largest to smallest eigenvalue}
#'     \item{rank}{Numerical rank of the matrix}
#'     \item{is_singular}{Logical: is the matrix numerically singular?}
#'     \item{warnings}{Character vector of any issues detected}
#'   }
#'
#' @details
#' For a proper maximum of the log-likelihood, the Hessian should be
#' negative definite (equivalently, -H should be positive definite).
#'
#' Common issues:
#' \itemize{
#'   \item **Not negative definite**: May indicate a saddle point rather than maximum
#'   \item **High condition number**: Indicates near-singularity, poorly identified parameters
#'   \item **Singular matrix**: Parameters are not identifiable from the data
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' x <- rnorm(100, mean = 5, sd = 2)
#'
#' result <- fit(
#'   function(mu, log_sigma) loglik_normal(mu, exp(log_sigma), x),
#'   params = c(mu = 0, log_sigma = 0)
#' )
#'
#' check <- check_hessian(result)
#' check
#'
#' # Access specific properties
#' check$is_negative_definite
#' check$condition_number
#' }
#'
#' @export
check_hessian <- function(object) {
  if (!inherits(object, "femtofit")) {
    stop("object must be a femtofit object")
  }

  H <- object$hessian
  if (is.null(H)) {
    stop("Hessian not available in fitted object")
  }

  # Observed information is -H
  I <- -H
  n <- nrow(I)

  # Compute eigenvalues
  eig <- eigen(I, symmetric = TRUE)
  eigenvalues <- eig$values

  # Check properties
  min_eig <- min(eigenvalues)
  max_eig <- max(eigenvalues)
  is_positive_definite <- all(eigenvalues > 0)

  # Condition number
  cond_num <- if (min_eig > 0) max_eig / min_eig else Inf

  # Numerical rank (eigenvalues above tolerance)
  tol <- max(dim(I)) * .Machine$double.eps * max_eig
  rank <- sum(eigenvalues > tol)

  is_singular <- rank < n

  # Collect warnings
  warnings <- character(0)

  if (!is_positive_definite) {
    warnings <- c(warnings,
      paste("Observed information is not positive definite.",
            "The Hessian has", sum(eigenvalues <= 0), "non-positive eigenvalue(s).",
            "This may indicate a saddle point rather than a maximum."))
  }

  if (cond_num > 1e10) {
    warnings <- c(warnings,
      paste("Very high condition number:", format(cond_num, scientific = TRUE),
            "Parameters may be poorly identified or nearly collinear."))
  } else if (cond_num > 1e6) {
    warnings <- c(warnings,
      paste("High condition number:", format(cond_num, scientific = TRUE),
            "Some parameters may be weakly identified."))
  }

  if (is_singular) {
    warnings <- c(warnings,
      paste("Matrix is numerically singular (rank", rank, "of", n, ").",
            "Some parameters may not be identifiable."))
  }

  structure(
    list(
      is_negative_definite = is_positive_definite,  # -H is pos def
      eigenvalues = eigenvalues,
      condition_number = cond_num,
      rank = rank,
      n_params = n,
      is_singular = is_singular,
      warnings = warnings
    ),
    class = "hessian_check"
  )
}


#' Check convergence diagnostics
#'
#' Examines convergence properties of a fitted model.
#'
#' @param object A `femtofit` object.
#' @param tol Tolerance for gradient norm (default 1e-4).
#'
#' @return A `convergence_check` object containing:
#'   \describe{
#'     \item{converged}{From the optimizer}
#'     \item{gradient_norm}{Norm of gradient at solution}
#'     \item{gradient_ok}{Is gradient norm below tolerance?}
#'     \item{iterations}{Number of iterations used}
#'     \item{loglik}{Final log-likelihood value}
#'     \item{warnings}{Character vector of any issues detected}
#'   }
#'
#' @details
#' At a maximum, the gradient should be (near) zero. Large gradient norms
#' may indicate premature termination or numerical issues.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' x <- rnorm(100, mean = 5, sd = 2)
#'
#' result <- fit(
#'   function(mu, log_sigma) loglik_normal(mu, exp(log_sigma), x),
#'   params = c(mu = 0, log_sigma = 0)
#' )
#'
#' check <- check_convergence(result)
#' check
#' }
#'
#' @export
check_convergence <- function(object, tol = 1e-4) {
  if (!inherits(object, "femtofit")) {
    stop("object must be a femtofit object")
  }

  g <- object$gradient
  if (is.null(g)) {
    g <- rep(NA, length(coef(object)))
  }

  # Compute gradient norm
  if (all(is.na(g))) {
    grad_norm <- NA
    gradient_ok <- NA
  } else {
    grad_norm <- sqrt(sum(g^2))
    gradient_ok <- grad_norm < tol
  }

  warnings <- character(0)

  if (!object$converged) {
    warnings <- c(warnings,
      "Optimizer did not report convergence. Results may be unreliable.")
  }

  if (!is.na(gradient_ok) && !gradient_ok) {
    warnings <- c(warnings,
      paste("Gradient norm", format(grad_norm, digits = 4),
            "exceeds tolerance", tol,
            ". Optimization may not have converged."))
  }

  if (!is.na(object$loglik) && !is.finite(object$loglik)) {
    warnings <- c(warnings,
      "Log-likelihood is not finite. Check for invalid parameter values.")
  }

  structure(
    list(
      converged = object$converged,
      gradient = g,
      gradient_norm = grad_norm,
      gradient_ok = gradient_ok,
      iterations = object$iterations,
      loglik = object$loglik,
      tolerance = tol,
      warnings = warnings
    ),
    class = "convergence_check"
  )
}


#' Comprehensive model diagnostics
#'
#' Runs all diagnostic checks on a fitted model.
#'
#' @param object A `femtofit` object.
#' @param ... Additional arguments passed to individual check functions.
#'
#' @return A `model_diagnostics` object containing results from all checks.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' x <- rnorm(100, mean = 5, sd = 2)
#'
#' result <- fit(
#'   function(mu, log_sigma) loglik_normal(mu, exp(log_sigma), x),
#'   params = c(mu = 0, log_sigma = 0)
#' )
#'
#' diagnostics(result)
#' }
#'
#' @export
diagnostics <- function(object, ...) {
  UseMethod("diagnostics")
}


#' @describeIn diagnostics Diagnostics for femtofit objects
#' @export
diagnostics.femtofit <- function(object, ...) {
  hessian_check <- check_hessian(object)
  convergence_check <- check_convergence(object, ...)

  # Collect all warnings
  all_warnings <- c(hessian_check$warnings, convergence_check$warnings)

  # Overall status
  overall_ok <- hessian_check$is_negative_definite &&
                !hessian_check$is_singular &&
                object$converged &&
                (is.na(convergence_check$gradient_ok) || convergence_check$gradient_ok)

  structure(
    list(
      hessian = hessian_check,
      convergence = convergence_check,
      ok = overall_ok,
      n_warnings = length(all_warnings),
      all_warnings = all_warnings
    ),
    class = "model_diagnostics"
  )
}


#' Print methods for diagnostic objects
#'
#' @param x A diagnostic object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.hessian_check <- function(x, ...) {
  cat("Hessian Diagnostics\n")
  cat("-------------------\n")
  cat("Negative definite:", x$is_negative_definite, "\n")
  cat("Condition number:", format(x$condition_number, scientific = TRUE), "\n")
  cat("Rank:", x$rank, "of", x$n_params, "\n")
  cat("Singular:", x$is_singular, "\n")

  if (length(x$warnings) > 0) {
    cat("\nWarnings:\n")
    for (w in x$warnings) {
      cat("  -", w, "\n")
    }
  } else {
    cat("\nNo issues detected.\n")
  }

  invisible(x)
}


#' @rdname print.hessian_check
#' @export
print.convergence_check <- function(x, ...) {
  cat("Convergence Diagnostics\n")
  cat("-----------------------\n")
  cat("Converged:", x$converged, "\n")
  cat("Iterations:", x$iterations, "\n")
  cat("Gradient norm:", format(x$gradient_norm, digits = 4), "\n")
  cat("Gradient OK:", x$gradient_ok, "(tolerance:", x$tolerance, ")\n")
  cat("Log-likelihood:", format(x$loglik, digits = 4), "\n")

  if (length(x$warnings) > 0) {
    cat("\nWarnings:\n")
    for (w in x$warnings) {
      cat("  -", w, "\n")
    }
  } else {
    cat("\nNo issues detected.\n")
  }

  invisible(x)
}


#' @rdname print.hessian_check
#' @export
print.model_diagnostics <- function(x, ...) {
  cat("Model Diagnostics\n")
  cat("=================\n\n")

  if (x$ok) {
    cat("Overall: OK\n\n")
  } else {
    cat("Overall: ISSUES DETECTED\n\n")
  }

  print(x$hessian)
  cat("\n")
  print(x$convergence)

  if (x$n_warnings > 0) {
    cat("\n-------------------\n")
    cat("Total warnings:", x$n_warnings, "\n")
  }

  invisible(x)
}


#' Check if standard errors are reliable
#'
#' A quick check of whether standard errors from the variance-covariance
#' matrix are likely to be reliable.
#'
#' @param object A `femtofit` object.
#'
#' @return Logical indicating whether SEs are likely reliable.
#'
#' @details
#' Standard errors may be unreliable if:
#' \itemize{
#'   \item The Hessian is not negative definite
#'   \item The model did not converge
#'   \item Any SE is NA or negative
#' }
#'
#' @export
se_reliable <- function(object) {
  if (!inherits(object, "femtofit")) {
    return(FALSE)
  }

  # Check vcov is valid
  V <- object$vcov
  if (is.null(V) || any(is.na(V))) {
    return(FALSE)
  }

  # Check SEs are positive and finite
  ses <- sqrt(diag(V))
  if (any(ses <= 0) || any(!is.finite(ses))) {
    return(FALSE)
  }

  # Check convergence
  if (!object$converged) {
    return(FALSE)
  }

  # Check Hessian is negative definite
  H <- object$hessian
  if (!is.null(H)) {
    eig <- eigen(-H, symmetric = TRUE, only.values = TRUE)$values
    if (any(eig <= 0)) {
      return(FALSE)
    }
  }

  TRUE
}
