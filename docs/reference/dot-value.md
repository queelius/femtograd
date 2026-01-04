# Subtraction for value objects

Element-wise subtraction or unary negation. Supports broadcasting when
one operand is a 1x1 scalar matrix.

## Usage

``` r
# S3 method for class 'value'
x - NULL
```

## Arguments

- x:

  A value object or numeric (minuend)

- y:

  A value object or numeric (subtrahend), or NULL for unary negation

## Value

A new value object representing x - y or -x
