# Generic function for the Backward pass for automatic differentiation (finds the gradient of every sub-node in the computational graph with respect to `e`). In other words, it is responsible for computing the gradient with respect to `e`.

This generic function should be implemented for specific classes. We
provide an implementation for `value` objects.

## Usage

``` r
backward(e, ...)
```

## Arguments

- e:

  An object for which the backward pass should be performed

- ...:

  additional arguments to pass
