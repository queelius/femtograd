# `value` object constructor

Creates a new value node in the computational graph with automatic
differentiation. Input is converted to matrix representation internally.

## Usage

``` r
val(data, label = "")
```

## Arguments

- data:

  Numeric value (scalar, vector, or matrix)

- label:

  Optional character label for debugging

## Value

A value object

## Details

All values are stored as matrices internally:

- `val(5)` creates a 1x1 matrix

- `val(c(1,2,3))` creates a 3x1 column vector

- `val(matrix(...))` preserves matrix dimensions
