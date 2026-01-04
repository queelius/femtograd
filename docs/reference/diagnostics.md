# Model Diagnostics

Functions for checking model fit quality, convergence, and potential
numerical issues.

Runs all diagnostic checks on a fitted model.

## Usage

``` r
diagnostics(object, ...)

# S3 method for class 'femtofit'
diagnostics(object, ...)
```

## Arguments

- object:

  A `femtofit` object.

- ...:

  Additional arguments passed to individual check functions.

## Value

A `model_diagnostics` object containing results from all checks.

## Methods (by class)

- `diagnostics(femtofit)`: Diagnostics for femtofit objects

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
x <- rnorm(100, mean = 5, sd = 2)

result <- fit(
  function(mu, log_sigma) loglik_normal(mu, exp(log_sigma), x),
  params = c(mu = 0, log_sigma = 0)
)

diagnostics(result)
} # }
```
