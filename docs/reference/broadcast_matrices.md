# Broadcast matrices for element-wise operations

Handles broadcasting when one operand is a 1x1 scalar matrix. Returns
the scalar value if broadcasting is needed, otherwise returns the
original matrix.

## Usage

``` r
broadcast_matrices(m1, m2)
```

## Arguments

- m1:

  First matrix

- m2:

  Second matrix

## Value

List with possibly broadcast matrices
