# Mean for value objects

Overrides base::mean to work with value objects. Computes the arithmetic
mean and properly propagates gradients (each input gets gradient / n).

## Usage

``` r
# S3 method for class 'value'
mean(x, ...)
```

## Arguments

- x:

  A value object or vector of value objects

- ...:

  Additional arguments (for compatibility with base::mean)

## Value

A new value object representing the mean
