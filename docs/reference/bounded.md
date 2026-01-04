# Transform to bounded interval

Maps unconstrained values to a bounded interval (lower, upper). Use this
for parameters with both lower and upper bounds.

## Usage

``` r
bounded(x, lower, upper)
```

## Arguments

- x:

  Unconstrained value (scalar, vector, or value object)

- lower:

  Lower bound of the interval

- upper:

  Upper bound of the interval

## Value

Value in (lower, upper)

## Details

The transformation is:

      bounded(x, a, b) = a + (b - a) * sigmoid(x)

This maps (-∞, ∞) to (a, b) smoothly.

## Examples

``` r
if (FALSE) { # \dontrun{
# Parameter in (0, 10)
raw <- val(0)
param <- bounded(raw, 0, 10)  # sigmoid(0) * 10 = 5
data(param)

# Correlation parameter in (-1, 1)
raw_rho <- val(1)
rho <- bounded(raw_rho, -1, 1)
} # }
```
