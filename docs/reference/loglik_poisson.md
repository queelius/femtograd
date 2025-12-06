# Poisson distribution log-likelihood

Computes the log-likelihood for i.i.d. Poisson observations. L(λ\|x) =
Σxᵢ*log(λ) - n*λ - Σlog(xᵢ!)

## Usage

``` r
loglik_poisson(lambda, x)
```

## Arguments

- lambda:

  Rate parameter λ (value object), must be positive

- x:

  Integer vector of observations (counts)

## Value

A value object representing the log-likelihood

## Details

The term Σlog(xᵢ!) is constant w.r.t. λ and is included for
completeness.

## Examples

``` r
if (FALSE) { # \dontrun{
x <- rpois(100, lambda = 3)
lambda <- val(1)
ll <- loglik_poisson(lambda, x)
backward(ll)
# MLE is mean(x)
} # }
```
