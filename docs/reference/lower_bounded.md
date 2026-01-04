# Transform to lower-bounded interval

Maps unconstrained values to (lower, ∞). Equivalent to lower +
positive(x).

## Usage

``` r
lower_bounded(x, lower)
```

## Arguments

- x:

  Unconstrained value

- lower:

  Lower bound

## Value

Value \> lower

## Details

The transformation is: lower_bounded(x, a) = a + exp(x)

## Examples

``` r
if (FALSE) { # \dontrun{
# Parameter > 2
raw <- val(0)
param <- lower_bounded(raw, 2)  # 2 + exp(0) = 3
} # }
```
