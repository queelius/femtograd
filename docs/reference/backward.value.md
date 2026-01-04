# Backward pass for value objects

Performs the backward pass (gradient computation) for a value object in
the computational graph with automatic differentiation. Initializes the
gradient of the output node to a matrix of ones (same shape as data).

## Usage

``` r
# S3 method for class 'value'
backward(x, ...)
```

## Arguments

- x:

  A value object for which the backward pass should be performed

- ...:

  pass additional arguments
