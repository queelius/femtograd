#' Optimization routines for maximum likelihood estimation
#'
#' These functions find the optimum (maximum or minimum) of differentiable
#' objective functions using femtograd's automatic differentiation.
#'
#' @name optimization
NULL


#' Gradient ascent/descent optimizer
#'
#' Finds the optimum of a function using gradient-based optimization.
#' Supports gradient clipping and adaptive step sizes.
#'
#' @param objective_fn Function taking list of value parameters, returns scalar
#' @param params List of value objects (initial parameter values)
#' @param lr Learning rate (step size), default 0.01
#' @param max_iter Maximum iterations, default 1000
#' @param tol Convergence tolerance on gradient norm, default 1e-6
#' @param maximize If TRUE (default), maximize; if FALSE, minimize
#' @param grad_clip Maximum gradient norm (NULL for no clipping)
#' @param verbose Print progress every N iterations (0 for silent)
#'
#' @return A list containing:
#'   \item{params}{List of value objects at optimum}
#'   \item{value}{Objective function value at optimum}
#'   \item{gradient}{Gradient at optimum}
#'   \item{iterations}{Number of iterations performed}
#'   \item{converged}{TRUE if gradient norm < tol}
#'
#' @examples
#' \dontrun{
#' # Find MLE for exponential distribution
#' x <- rexp(100, rate = 2)
#' loglik <- function(p) loglik_exponential(p[[1]], x)
#' result <- gradient_ascent(loglik, list(val(1)))
#' data(result$params[[1]])  # Should be close to 1/mean(x)
#' }
#'
#' @export
gradient_ascent <- function(objective_fn, params, lr = 0.01, max_iter = 1000,
                            tol = 1e-6, maximize = TRUE, grad_clip = NULL,
                            verbose = 0)
{
  n <- length(params)
  sign <- if (maximize) 1 else -1

  # Extract initial values
  param_values <- sapply(params, data)

  converged <- FALSE
  iter <- 0
  obj_value <- NA
  grad_vec <- rep(NA, n)

  for (iter in seq_len(max_iter)) {
    # Create fresh value objects
    current_params <- lapply(param_values, val)

    # Compute objective and gradients
    obj <- objective_fn(current_params)
    obj_value <- data(obj)
    backward(obj)

    # Extract gradients
    grad_vec <- sapply(current_params, grad)
    grad_norm <- sqrt(sum(grad_vec^2))

    # Check convergence
    if (grad_norm < tol) {
      converged <- TRUE
      break
    }

    # Gradient clipping
    if (!is.null(grad_clip) && grad_norm > grad_clip) {
      grad_vec <- grad_vec * (grad_clip / grad_norm)
    }

    # Update parameters
    param_values <- param_values + sign * lr * grad_vec

    # Progress report
    if (verbose > 0 && iter %% verbose == 0) {
      cat(sprintf("Iter %d: objective = %.6f, grad_norm = %.6e\n",
                  iter, obj_value, grad_norm))
    }
  }

  # Return final params as value objects
  final_params <- lapply(param_values, val)

  list(
    params = final_params,
    value = obj_value,
    gradient = grad_vec,
    iterations = iter,
    converged = converged
  )
}


#' Gradient descent (minimize)
#'
#' Convenience wrapper for gradient_ascent with maximize=FALSE.
#'
#' @inheritParams gradient_ascent
#' @export
gradient_descent <- function(objective_fn, params, lr = 0.01, max_iter = 1000,
                             tol = 1e-6, grad_clip = NULL, verbose = 0)
{
  gradient_ascent(objective_fn, params, lr = lr, max_iter = max_iter,
                  tol = tol, maximize = FALSE, grad_clip = grad_clip,
                  verbose = verbose)
}


#' Newton-Raphson optimizer
#'
#' Finds the optimum using Newton-Raphson method with exact Hessian.
#' Uses second-order information for faster convergence near optimum.
#'
#' @param objective_fn Function taking list of value parameters, returns scalar
#' @param params List of value objects (initial parameter values)
#' @param max_iter Maximum iterations, default 100
#' @param tol Convergence tolerance on step size, default 1e-8
#' @param maximize If TRUE (default), maximize; if FALSE, minimize
#' @param step_scale Scale factor for Newton step (< 1 for damping), default 1
#' @param verbose Print progress every N iterations (0 for silent)
#'
#' @return A list containing:
#'   \item{params}{List of value objects at optimum}
#'   \item{value}{Objective function value at optimum}
#'   \item{gradient}{Gradient at optimum}
#'   \item{hessian}{Hessian at optimum}
#'   \item{iterations}{Number of iterations performed}
#'   \item{converged}{TRUE if step size < tol}
#'
#' @details
#' Newton-Raphson update: θ_{n+1} = θ_n - H⁻¹ g
#' For maximization, uses: θ_{n+1} = θ_n - H⁻¹ g (H is negative definite)
#' For minimization, uses: θ_{n+1} = θ_n - H⁻¹ g (H is positive definite)
#'
#' The Hessian is computed via forward-over-reverse AD at each iteration.
#' This is exact but can be slow for many parameters.
#'
#' @examples
#' \dontrun{
#' # Find MLE for normal distribution
#' x <- rnorm(100, mean = 5, sd = 2)
#' loglik <- function(p) loglik_normal(p[[1]], p[[2]], x)
#' result <- newton_raphson(loglik, list(val(0), val(1)))
#' sapply(result$params, data)  # Should be close to c(mean(x), sd(x))
#' }
#'
#' @export
newton_raphson <- function(objective_fn, params, max_iter = 100, tol = 1e-8,
                           maximize = TRUE, step_scale = 1, verbose = 0)
{
  n <- length(params)
  sign <- if (maximize) 1 else -1

  # Extract initial values
  param_values <- sapply(params, data)

  converged <- FALSE
  iter <- 0
  obj_value <- NA
  grad_vec <- rep(NA, n)
  H <- matrix(NA, n, n)

  for (iter in seq_len(max_iter)) {
    # Create fresh value objects
    current_params <- lapply(param_values, val)

    # Compute objective, gradient, and Hessian
    obj_value <- data(objective_fn(current_params))
    grad_vec <- gradient(objective_fn, current_params)
    H <- hessian(objective_fn, current_params)

    # Newton step: -H⁻¹ g
    # For maximization of log-likelihood, H is typically negative definite
    # So we use: θ_{n+1} = θ_n - H⁻¹ g
    tryCatch({
      step <- solve(H, grad_vec)
    }, error = function(e) {
      warning("Hessian is singular, using gradient step instead")
      step <- sign * 0.01 * grad_vec
    })

    step <- -step * step_scale

    # Check convergence
    step_norm <- sqrt(sum(step^2))
    if (step_norm < tol) {
      converged <- TRUE
      break
    }

    # Update parameters
    param_values <- param_values + step

    # Progress report
    if (verbose > 0 && iter %% verbose == 0) {
      cat(sprintf("Iter %d: objective = %.6f, step_norm = %.6e\n",
                  iter, obj_value, step_norm))
    }
  }

  # Return final params as value objects
  final_params <- lapply(param_values, val)

  list(
    params = final_params,
    value = obj_value,
    gradient = grad_vec,
    hessian = H,
    iterations = iter,
    converged = converged
  )
}


#' Fisher scoring optimizer
#'
#' Similar to Newton-Raphson but uses expected Fisher information
#' (negative expected Hessian) instead of observed Hessian.
#' More stable for some problems.
#'
#' For regular exponential families, Fisher scoring is equivalent to
#' Newton-Raphson since observed = expected information.
#'
#' @param loglik_fn Log-likelihood function
#' @param params List of value objects (initial parameter values)
#' @param max_iter Maximum iterations, default 100
#' @param tol Convergence tolerance, default 1e-8
#' @param verbose Print progress every N iterations (0 for silent)
#'
#' @return Same structure as newton_raphson
#'
#' @details
#' Fisher scoring: θ_{n+1} = θ_n + I(θ_n)⁻¹ S(θ_n)
#' where I = -E[H] (Fisher information) and S = gradient (score).
#'
#' This implementation uses the observed Hessian as an approximation
#' to the expected Hessian, making it identical to Newton-Raphson.
#' For a true Fisher scoring implementation, one would need to compute
#' E[H] analytically or via Monte Carlo.
#'
#' @export
fisher_scoring <- function(loglik_fn, params, max_iter = 100, tol = 1e-8,
                           verbose = 0)
{
  newton_raphson(loglik_fn, params, max_iter = max_iter, tol = tol,
                 maximize = TRUE, verbose = verbose)
}


#' Find MLE with standard errors
#'
#' Convenience function that finds the MLE and computes standard errors
#' from the observed Fisher information.
#'
#' @param loglik_fn Log-likelihood function taking list of value parameters
#' @param params List of value objects (initial parameter values)
#' @param method Optimization method: "newton" or "gradient"
#' @param ... Additional arguments passed to optimizer
#'
#' @return A list containing:
#'   \item{estimate}{MLE as numeric vector}
#'   \item{se}{Standard errors}
#'   \item{vcov}{Variance-covariance matrix}
#'   \item{loglik}{Log-likelihood at MLE}
#'   \item{hessian}{Hessian at MLE}
#'   \item{converged}{Convergence indicator}
#'
#' @examples
#' \dontrun{
#' x <- rnorm(100, mean = 5, sd = 2)
#' loglik <- function(p) loglik_normal(p[[1]], p[[2]], x)
#' result <- find_mle(loglik, list(val(0), val(1)))
#' result$estimate  # MLEs
#' result$se        # Standard errors
#' }
#'
#' @export
find_mle <- function(loglik_fn, params, method = "newton", ...)
{
  # Run optimizer
  if (method == "newton") {
    opt <- newton_raphson(loglik_fn, params, maximize = TRUE, ...)
  } else if (method == "gradient") {
    opt <- gradient_ascent(loglik_fn, params, maximize = TRUE, ...)
    # Need to compute Hessian for standard errors
    opt$hessian <- hessian(loglik_fn, opt$params)
  } else {
    stop("method must be 'newton' or 'gradient'")
  }

  # Extract estimates
  estimate <- sapply(opt$params, data)

  # Compute standard errors from Hessian
  vcov <- tryCatch({
    solve(-opt$hessian)
  }, error = function(e) {
    warning("Could not invert Hessian for standard errors")
    matrix(NA, length(params), length(params))
  })

  se <- sqrt(diag(vcov))

  list(
    estimate = estimate,
    se = se,
    vcov = vcov,
    loglik = opt$value,
    hessian = opt$hessian,
    gradient = opt$gradient,
    converged = opt$converged,
    iterations = opt$iterations
  )
}


#' Compute confidence intervals from MLE results
#'
#' @param mle_result Result from find_mle
#' @param level Confidence level, default 0.95
#' @param param_names Optional names for parameters
#'
#' @return A matrix with columns: estimate, se, lower, upper
#'
#' @export
confint_mle <- function(mle_result, level = 0.95, param_names = NULL)
{
  z <- qnorm(1 - (1 - level) / 2)
  lower <- mle_result$estimate - z * mle_result$se
  upper <- mle_result$estimate + z * mle_result$se

  result <- cbind(
    estimate = mle_result$estimate,
    se = mle_result$se,
    lower = lower,
    upper = upper
  )

  if (!is.null(param_names)) {
    rownames(result) <- param_names
  }

  result
}


#' Wald test for hypothesis testing
#'
#' Computes Wald test statistic for H0: θ = θ₀
#' W = (θ̂ - θ₀)' I(θ̂) (θ̂ - θ₀) ~ χ²(k)
#'
#' @param mle_result Result from find_mle
#' @param null_values Null hypothesis values (default: 0 for all params)
#'
#' @return A list with test statistic, df, and p-value
#'
#' @export
wald_test <- function(mle_result, null_values = NULL)
{
  k <- length(mle_result$estimate)
  if (is.null(null_values)) null_values <- rep(0, k)

  diff <- mle_result$estimate - null_values

  # Wald statistic: (θ̂ - θ₀)' Σ⁻¹ (θ̂ - θ₀)
  # where Σ = Var(θ̂) = -H⁻¹
  W <- tryCatch({
    t(diff) %*% solve(mle_result$vcov) %*% diff
  }, error = function(e) {
    NA
  })

  list(
    statistic = as.numeric(W),
    df = k,
    p_value = 1 - pchisq(as.numeric(W), df = k)
  )
}
