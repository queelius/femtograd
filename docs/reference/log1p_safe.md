# Log1p with underflow protection

Computes log(1 + x) with protection for very small x.

## Usage

``` r
log1p_safe(x)
```

## Arguments

- x:

  A value object or numeric

## Value

log(1 + x), accurate even for x near zero

## Details

R's built-in log1p is already stable, but this version adds autodiff
support and handles the case where 1 + x might underflow to 1 for very
small x.
