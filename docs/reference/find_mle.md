# Find MLE with standard errors

Convenience function that finds the MLE and computes standard errors
from the observed Fisher information.

## Usage

``` r
find_mle(loglik_fn, params, method = "newton", ...)
```

## Arguments

- loglik_fn:

  Log-likelihood function taking list of value parameters

- params:

  List of value objects (initial parameter values)

- method:

  Optimization method: "newton" or "gradient"

- ...:

  Additional arguments passed to optimizer

## Value

A list containing:

- estimate:

  MLE as numeric vector

- se:

  Standard errors

- vcov:

  Variance-covariance matrix

- loglik:

  Log-likelihood at MLE

- hessian:

  Hessian at MLE

- converged:

  Convergence indicator

## Examples

``` r
if (FALSE) { # \dontrun{
x <- rnorm(100, mean = 5, sd = 2)
loglik <- function(p) loglik_normal(p[[1]], p[[2]], x)
result <- find_mle(loglik, list(val(0), val(1)))
result$estimate  # MLEs
result$se        # Standard errors
} # }
```
