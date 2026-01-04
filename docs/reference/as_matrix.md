# Convert input to matrix representation

Internal helper to ensure all data is stored as matrices.

- Scalars become 1x1 matrices

- Vectors become n x 1 column matrices

- Matrices remain unchanged

## Usage

``` r
as_matrix(x)
```

## Arguments

- x:

  Numeric scalar, vector, or matrix

## Value

A matrix
