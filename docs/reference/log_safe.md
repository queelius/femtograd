# Safe logarithm (handles zeros)

Computes log(x) with protection against log(0) = -Inf. Optionally clamps
input to a minimum value.

## Usage

``` r
log_safe(x, eps = .Machine$double.eps)
```

## Arguments

- x:

  A value object or numeric

- eps:

  Minimum value to clamp x to (default: .Machine\$double.eps)

## Value

log(max(x, eps))

## Details

This is useful when computing log-likelihoods where probabilities might
numerically become zero. The gradient is computed as if the clamping
didn't happen (straight-through estimator) when x \> eps.
