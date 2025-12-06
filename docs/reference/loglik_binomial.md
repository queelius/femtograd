# Binomial distribution log-likelihood

Computes the log-likelihood for binomial observations. L(p\|x,n) =
Σxᵢ\*log(p) + (nᵢ-xᵢ)\*log(1-p) + log(C(nᵢ,xᵢ))

## Usage

``` r
loglik_binomial(p, x, size)
```

## Arguments

- p:

  Success probability (value object), must be in (0,1)

- x:

  Integer vector of successes

- size:

  Integer vector of trial counts (or single value if constant)

## Value

A value object representing the log-likelihood

## Examples

``` r
if (FALSE) { # \dontrun{
x <- rbinom(50, size = 10, prob = 0.3)
p <- val(0.5)
ll <- loglik_binomial(p, x, size = 10)
backward(ll)
# MLE is sum(x) / sum(size)
} # }
```
