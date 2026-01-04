# Addition for value objects

Element-wise addition of two value objects or a value and numeric. Both
operands are converted to matrices internally. Supports broadcasting
when one operand is a 1x1 scalar matrix.

## Usage

``` r
# S3 method for class 'value'
e1 + e2
```

## Arguments

- e1:

  A value object or numeric

- e2:

  A value object or numeric

## Value

A new value object representing the element-wise sum
