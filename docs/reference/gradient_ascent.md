# Gradient ascent/descent optimizer

Finds the optimum of a function using gradient-based optimization.
Supports gradient clipping and adaptive step sizes.

## Usage

``` r
gradient_ascent(
  objective_fn,
  params,
  lr = 0.01,
  max_iter = 1000,
  tol = 1e-06,
  maximize = TRUE,
  grad_clip = NULL,
  verbose = 0
)
```

## Arguments

- objective_fn:

  Function taking list of value parameters, returns scalar

- params:

  List of value objects (initial parameter values)

- lr:

  Learning rate (step size), default 0.01

- max_iter:

  Maximum iterations, default 1000

- tol:

  Convergence tolerance on gradient norm, default 1e-6

- maximize:

  If TRUE (default), maximize; if FALSE, minimize

- grad_clip:

  Maximum gradient norm (NULL for no clipping)

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

- iterations:

  Number of iterations performed

- converged:

  TRUE if gradient norm \< tol

## Examples

``` r
if (FALSE) { # \dontrun{
# Find MLE for exponential distribution
x <- rexp(100, rate = 2)
loglik <- function(p) loglik_exponential(p[[1]], x)
result <- gradient_ascent(loglik, list(val(1)))
data(result$params[[1]])  # Should be close to 1/mean(x)
} # }
```
