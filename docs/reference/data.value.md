# Retrieve the data matrix from a value object

Retrieve the data matrix from a value object

## Usage

``` r
# S3 method for class 'value'
data(x, drop = TRUE, ...)
```

## Arguments

- x:

  A value object

- drop:

  If TRUE (default) and result is 1x1, return scalar. Set to FALSE to
  always return a matrix.

- ...:

  additional arguments to pass

## Value

The data as scalar (if 1x1 and drop=TRUE) or matrix
