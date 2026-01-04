# BFGS quasi-Newton optimizer

Finds the optimum using the BFGS algorithm, which approximates the
inverse Hessian using gradient information. More efficient than
Newton-Raphson when Hessian computation is expensive.

## Usage

``` r
bfgs(
  objective_fn,
  params,
  max_iter = 1000,
  tol = 1e-06,
  maximize = FALSE,
  line_search_fn = NULL,
  verbose = 0
)
```

## Arguments

- objective_fn:

  Function taking list of value parameters, returns scalar

- params:

  List of value objects (initial parameter values)

- max_iter:

  Maximum iterations, default 1000

- tol:

  Convergence tolerance on gradient norm, default 1e-6

- maximize:

  If TRUE, maximize; if FALSE (default), minimize

- line_search_fn:

  Line search function (default: backtracking)

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

- inv_hessian:

  Approximate inverse Hessian

- iterations:

  Number of iterations performed

- converged:

  TRUE if gradient norm \< tol

## Details

BFGS maintains an approximation B to the inverse Hessian, updated as:
B_k+1 = (I - ρ*s*y') B_k (I - ρ*y*s') + ρ*s*s' where s = x_k+1 - x_k, y
= g_k+1 - g_k, ρ = 1/(y'\*s)

The search direction is d = -B\*g, followed by line search.

## Examples

``` r
if (FALSE) { # \dontrun{
# Minimize Rosenbrock function
rosenbrock <- function(p) {
  x <- p[[1]]; y <- p[[2]]
  (1 - x)^2 + 100 * (y - x^2)^2
}
result <- bfgs(rosenbrock, list(val(-1), val(-1)))
sapply(result$params, data)  # Should be close to c(1, 1)
} # }
```
