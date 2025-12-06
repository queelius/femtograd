# Summation for value objects

Overrides base::sum to work with value objects. Supports mixed lists of
value objects and numeric scalars.

## Usage

``` r
# S3 method for class 'value'
sum(..., na.rm = FALSE)
```

## Arguments

- ...:

  value objects and/or numeric scalars to sum

- na.rm:

  Logical, whether to remove NA values (passed to base sum for scalars)

## Value

A new value object representing the sum
