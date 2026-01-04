# Weibull distribution log-likelihood

Computes the log-likelihood for i.i.d. Weibull observations. L(k, λ \|
x) = n*log(k) - n*k\*log(λ) + (k-1)\*Σlog(xᵢ) - Σ(xᵢ/λ)^k

## Usage

``` r
loglik_weibull(shape, scale, x)
```

## Arguments

- shape:

  Shape parameter k (value object), must be positive

- scale:

  Scale parameter λ (value object), must be positive

- x:

  Numeric vector of observations (must be positive)

## Value

A value object representing the log-likelihood

## Details

The Weibull distribution is commonly used in survival analysis and
reliability engineering. It generalizes the exponential distribution
(k=1 gives exponential with rate 1/λ).

## Examples

``` r
if (FALSE) { # \dontrun{
x <- rweibull(100, shape = 2, scale = 3)

# Using log-parameterization for positivity
result <- fit(
  function(log_shape, log_scale) {
    shape <- exp(log_shape)
    scale <- exp(log_scale)
    loglik_weibull(shape, scale, x)
  },
  params = c(log_shape = 0, log_scale = 0)
)
} # }
```
