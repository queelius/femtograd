# Bernoulli distribution log-likelihood

Special case of binomial with size=1. L(p\|x) = Σxᵢ\*log(p) +
(1-xᵢ)\*log(1-p)

## Usage

``` r
loglik_bernoulli(p, x)
```

## Arguments

- p:

  Success probability (value object), must be in (0,1)

- x:

  Binary vector (0 or 1)

## Value

A value object representing the log-likelihood
