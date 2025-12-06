# Fisher scoring optimizer

Similar to Newton-Raphson but uses expected Fisher information (negative
expected Hessian) instead of observed Hessian. More stable for some
problems.

## Usage

``` r
fisher_scoring(loglik_fn, params, max_iter = 100, tol = 1e-08, verbose = 0)
```

## Arguments

- loglik_fn:

  Log-likelihood function

- params:

  List of value objects (initial parameter values)

- max_iter:

  Maximum iterations, default 100

- tol:

  Convergence tolerance, default 1e-8

- verbose:

  Print progress every N iterations (0 for silent)

## Value

Same structure as newton_raphson

## Details

For regular exponential families, Fisher scoring is equivalent to
Newton-Raphson since observed = expected information.

Fisher scoring: θ_n+1 = θ_n + I(θ_n)⁻¹ S(θ_n) where I = -EH (Fisher
information) and S = gradient (score).

This implementation uses the observed Hessian as an approximation to the
expected Hessian, making it identical to Newton-Raphson. For a true
Fisher scoring implementation, one would need to compute EH analytically
or via Monte Carlo.
