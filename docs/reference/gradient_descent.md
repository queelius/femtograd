# Gradient descent (minimize)

Convenience wrapper for gradient_ascent with maximize=FALSE.

## Usage

``` r
gradient_descent(
  objective_fn,
  params,
  lr = 0.01,
  max_iter = 1000,
  tol = 1e-06,
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

- grad_clip:

  Maximum gradient norm (NULL for no clipping)

- verbose:

  Print progress every N iterations (0 for silent)
