# Newton-Raphson optimizer

Finds the optimum using Newton-Raphson method with exact Hessian. Uses
second-order information for faster convergence near optimum.

## Usage

``` r
newton_raphson(
  objective_fn,
  params,
  max_iter = 100,
  tol = 1e-08,
  maximize = TRUE,
  step_scale = 1,
  verbose = 0
)
```

## Arguments

- objective_fn:

  Function taking list of value parameters, returns scalar

- params:

  List of value objects (initial parameter values)

- max_iter:

  Maximum iterations, default 100

- tol:

  Convergence tolerance on step size, default 1e-8

- maximize:

  If TRUE (default), maximize; if FALSE, minimize

- step_scale:

  Scale factor for Newton step (\< 1 for damping), default 1

- verbose:

  Print progress every N iterations (0 for silent)

## Value

A list containing:

- params:

  List of value objects at optimum

- value:

  Objective function value at optimum

- gradient:

  Gradient at optimum

- hessian:

  Hessian at optimum

- iterations:

  Number of iterations performed

- converged:

  TRUE if step size \< tol

## Details

Newton-Raphson update: θ_n+1 = θ_n - H⁻¹ g For maximization, uses: θ_n+1
= θ_n - H⁻¹ g (H is negative definite) For minimization, uses: θ_n+1 =
θ_n - H⁻¹ g (H is positive definite)

The Hessian is computed via forward-over-reverse AD at each iteration.
This is exact but can be slow for many parameters.

## Examples

``` r
if (FALSE) { # \dontrun{
# Find MLE for normal distribution
x <- rnorm(100, mean = 5, sd = 2)
loglik <- function(p) loglik_normal(p[[1]], p[[2]], x)
result <- newton_raphson(loglik, list(val(0), val(1)))
sapply(result$params, data)  # Should be close to c(mean(x), sd(x))
} # }
```
