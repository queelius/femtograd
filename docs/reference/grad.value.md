# Gradient of a `value` object `x` with respect to `e` in `backward(e)`, e.g., dx/de. (applies the chain rule)

Gradient of a `value` object `x` with respect to `e` in `backward(e)`,
e.g., dx/de. (applies the chain rule)

## Usage

``` r
# S3 method for class 'value'
grad(x, drop = TRUE, ...)
```

## Arguments

- x:

  A value object

- drop:

  If TRUE (default) and result is 1x1, return scalar. Set to FALSE to
  always return a matrix.

- ...:

  pass additional arguments

## Value

The gradient as scalar (if 1x1 and drop=TRUE) or matrix
