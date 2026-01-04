# Stable sigmoid function

Computes sigmoid(x) = 1/(1+exp(-x)) with overflow protection.

## Usage

``` r
sigmoid_stable(x)
```

## Arguments

- x:

  A value object or numeric

## Value

Sigmoid values in (0, 1)

## Details

For large positive x, exp(-x) underflows to 0, giving sigmoid = 1
(correct). For large negative x, we use sigmoid(x) = exp(x)/(1+exp(x))
to avoid overflow.
