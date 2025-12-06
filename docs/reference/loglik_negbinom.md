# Negative binomial log-likelihood

Computes the log-likelihood for negative binomial (failures until r
successes). L(r,p\|x) = Σ

## Usage

``` r
loglik_negbinom(r, p, x)
```

## Arguments

- r:

  Number of successes parameter (value object or fixed positive)

- p:

  Success probability (value object), must be in (0,1)

- x:

  Integer vector of failure counts

## Value

A value object representing the log-likelihood
