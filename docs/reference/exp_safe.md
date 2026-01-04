# Stable exp function (with overflow protection)

Computes exp(x) with optional clamping to prevent Inf.

## Usage

``` r
exp_safe(x, max_val = 709)
```

## Arguments

- x:

  A value object or numeric

- max_val:

  Maximum value for x to prevent overflow (default: 709)

## Value

exp(min(x, max_val))

## Details

In double precision, exp(710) overflows to Inf. This function clamps the
input to prevent overflow while maintaining correct gradients.
