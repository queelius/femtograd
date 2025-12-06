# Backward pass for value objects

Performs the backward pass (gradient computation) for a value object in
the computational graph with automatic differentiation.

## Usage

``` r
# S3 method for class 'value'
backward(x)
```

## Arguments

- e:

  A value object for which the backward pass should be performed

- ...:

  pass additional arguments
