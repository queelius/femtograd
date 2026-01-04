# L-BFGS optimizer (limited memory BFGS)

Memory-efficient variant of BFGS that stores only the last m correction
pairs instead of the full inverse Hessian approximation. Suitable for
large-scale optimization.

## Usage

``` r
lbfgs(
  objective_fn,
  params,
  m = 10,
  max_iter = 1000,
  tol = 1e-06,
  maximize = FALSE,
  verbose = 0
)
```

## Arguments

- objective_fn:

  Function taking list of value parameters, returns scalar

- params:

  List of value objects (initial parameter values)

- m:

  Number of correction pairs to store, default 10

- max_iter:

  Maximum iterations, default 1000

- tol:

  Convergence tolerance on gradient norm, default 1e-6

- maximize:

  If TRUE, maximize; if FALSE (default), minimize

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

## Details

L-BFGS uses the two-loop recursion algorithm to compute the search
direction without explicitly forming the inverse Hessian. Memory usage
is O(m\*n) instead of O(n^2) for full BFGS.
