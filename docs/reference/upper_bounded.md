# Transform to upper-bounded interval

Maps unconstrained values to (-∞, upper). Equivalent to upper -
positive(x).

## Usage

``` r
upper_bounded(x, upper)
```

## Arguments

- x:

  Unconstrained value

- upper:

  Upper bound

## Value

Value \< upper

## Details

The transformation is: upper_bounded(x, b) = b - exp(x)

## Examples

``` r
if (FALSE) { # \dontrun{
# Parameter < 10
raw <- val(0)
param <- upper_bounded(raw, 10)  # 10 - exp(0) = 9
} # }
```
