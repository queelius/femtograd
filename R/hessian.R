#' Compute Hessian matrix via forward-over-reverse automatic differentiation
#'
#' Computes the Hessian matrix (matrix of second partial derivatives) of a
#' scalar function with respect to a vector of parameters. Uses forward-mode
#' AD on top of reverse-mode AD for efficient, accurate computation.
#'
#' @param loss_fn A function that takes a list of parameters and returns a
#'   scalar loss/objective. The function must use femtograd operations.
#' @param params A list of value objects (the parameters)
#' @param value_creator Function to create value objects (default: val).
#'   This allows customization for different value types.
#'
#' @return A numeric matrix of dimension (n x n) where n = length(params),
#'   containing the Hessian d²f/dθᵢdθⱼ
#'
#' @details
#' The Hessian is computed using forward-over-reverse mode AD:
#' \enumerate{
#'   \item For each parameter i, create dual numbers where:
#'     \itemize{
#'       \item Primal = value object holding parameter value
#'       \item Tangent = value object (1 for param i, 0 for others)
#'     }
#'   \item Run the loss function with these dual-value inputs
#'   \item The tangent of the loss is df/dθᵢ (as a value expression)
#'   \item Call backward() on the tangent loss
#'   \item Each tangent parameter's gradient gives d²f/dθⱼdθᵢ
#' }
#'
#' This is the "forward-over-reverse" pattern: forward-mode (dual numbers)
#' computes df/dθᵢ symbolically, then reverse-mode differentiates that
#' to get d²f/dθⱼdθᵢ.
#'
#' @examples
#' \dontrun{
#' # Hessian of f(x,y) = x^2 + x*y + y^2
#' loss_fn <- function(p) {
#'   x <- p[[1]]
#'   y <- p[[2]]
#'   x^2 + x*y + y^2
#' }
#' params <- list(val(1), val(2))
#' H <- hessian(loss_fn, params)
#' # H should be [[2, 1], [1, 2]]
#' }
#'
#' @export
hessian <- function(loss_fn, params, value_creator = val)
{
  n <- length(params)
  H <- matrix(0, nrow = n, ncol = n)

  # Extract current parameter values (as scalars using drop=TRUE)
  param_values <- sapply(params, function(p) data(p, drop = TRUE))

  for (i in seq_len(n)) {
    # Create fresh value objects with dual structure
    # Primal: value objects for the function value
    # Tangent: value objects where param i has tangent 1, others 0
    # This computes df/dθᵢ via forward-mode as a value expression
    dual_params <- lapply(seq_len(n), function(j) {
      primal_val <- value_creator(param_values[j])
      tangent_val <- value_creator(if (j == i) 1 else 0)
      dual$new(primal_val, tangent_val)
    })

    # Forward pass: compute loss with dual-value inputs
    # primal(dual_loss) = f(θ) as a value expression
    # tangent(dual_loss) = df/dθᵢ as a value expression
    dual_loss <- loss_fn(dual_params)

    # Extract tangent loss: this is df/dθᵢ represented as a value graph
    tangent_loss <- tangent(dual_loss)

    # Run backward on tangent loss to get d/dθⱼ(df/dθᵢ) = d²f/dθⱼdθᵢ
    # This is reverse-mode on the tangent expression
    backward(tangent_loss)

    # The PRIMAL parameter's gradient gives the Hessian entry
    # (tangent expression depends on both primal and tangent params;
    #  differentiating w.r.t. primal gives the second derivative)
    # Note: When a scalar parameter is broadcast across vector operations,
    # grad() returns a matrix - sum all elements to get the total derivative.
    for (j in seq_len(n)) {
      primal_param <- primal(dual_params[[j]])
      g <- grad(primal_param)
      H[j, i] <- sum(g)  # Sum all elements of the gradient matrix
    }
  }

  H
}


#' Compute gradient as a numeric vector
#'
#' Convenience function to compute the gradient of a loss function
#' at the current parameter values.
#'
#' @param loss_fn A function taking a list of value parameters
#' @param params A list of value objects
#'
#' @return A numeric vector of gradients (one per parameter)
#' @export
gradient <- function(loss_fn, params)
{
  # Create fresh value objects (extract scalar values with drop=TRUE)
  param_values <- sapply(params, function(p) data(p, drop = TRUE))
  fresh_params <- lapply(param_values, val)

  # Compute loss and backward
  loss <- loss_fn(fresh_params)
  backward(loss)

  # Extract gradients (sum each gradient matrix to get scalar)
  sapply(fresh_params, function(p) sum(grad(p)))
}


#' Compute observed Fisher information matrix
#'
#' For a log-likelihood function, the observed Fisher information is
#' the negative Hessian evaluated at the MLE. This is used for computing
#' standard errors of MLEs.
#'
#' @param loglik_fn Log-likelihood function taking a list of value parameters
#' @param params List of value objects at MLE
#'
#' @return The observed Fisher information matrix (negative Hessian)
#'
#' @details
#' The observed Fisher information I(θ̂) = -H(θ̂) where H is the Hessian
#' of the log-likelihood. Standard errors are sqrt(diag(I⁻¹)).
#'
#' @export
fisher_information <- function(loglik_fn, params)
{
  -hessian(loglik_fn, params)
}


#' Compute standard errors from Hessian
#'
#' Extracts standard errors from the Hessian of a log-likelihood.
#' SE(θ̂ᵢ) = sqrt([I(θ̂)⁻¹]ᵢᵢ) = sqrt([-H⁻¹]ᵢᵢ)
#'
#' @param hess The Hessian matrix (from \code{hessian()})
#' @param is_loglik If TRUE, treats hess as Hessian of log-likelihood
#'   and uses -H⁻¹. If FALSE, uses H⁻¹ directly.
#'
#' @return A numeric vector of standard errors
#' @export
std_errors <- function(hess, is_loglik = TRUE)
{
  if (is_loglik) {
    # For log-likelihood, variance is -H⁻¹
    vcov <- solve(-hess)
  } else {
    vcov <- solve(hess)
  }
  sqrt(diag(vcov))
}


#' Compute variance-covariance matrix from Hessian
#'
#' @param hess The Hessian matrix
#' @param is_loglik If TRUE, returns -H⁻¹ (for log-likelihood)
#'
#' @return The variance-covariance matrix
#' @export
vcov_matrix <- function(hess, is_loglik = TRUE)
{
  if (is_loglik) {
    solve(-hess)
  } else {
    solve(hess)
  }
}
