# Exponential distribution log-likelihood

Computes the log-likelihood for i.i.d. exponential observations. L(λ\|x)
= n*log(λ) - λ*Σxᵢ

## Usage

``` r
loglik_exponential(rate, x)
```

## Arguments

- rate:

  Rate parameter λ (value object), must be positive

- x:

  Numeric vector of observations (must be non-negative)

## Value

A value object representing the log-likelihood

## Examples

``` r
if (FALSE) { # \dontrun{
x <- rexp(50, rate = 2)
rate <- val(1)
ll <- loglik_exponential(rate, x)
backward(ll)
grad(rate)  # should be n/rate - sum(x)
} # }
```
