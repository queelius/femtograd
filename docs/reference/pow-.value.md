# Power operation for value objects.

Element-wise power operation. Supports broadcasting when one operand is
a 1x1 scalar matrix.

## Usage

``` r
# S3 method for class 'value'
b^e
```

## Arguments

- b:

  A value object or numeric (base)

- e:

  A value object or numeric (exponent)

## Value

A new value object representing the element-wise power b^e
