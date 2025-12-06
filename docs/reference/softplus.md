# Softplus function for value objects

Computes log(1 + exp(x)), a smooth approximation to max(0, x). Useful
for ensuring positivity of parameters (e.g., variances).

## Usage

``` r
softplus(x)
```

## Arguments

- x:

  A value object

## Value

A new value object representing softplus(x)
