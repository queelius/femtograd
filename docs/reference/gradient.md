# Compute gradient as a numeric vector

Convenience function to compute the gradient of a loss function at the
current parameter values.

## Usage

``` r
gradient(loss_fn, params)
```

## Arguments

- loss_fn:

  A function taking a list of value parameters

- params:

  A list of value objects

## Value

A numeric vector of gradients
