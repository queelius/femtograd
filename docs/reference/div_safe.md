# Safe division (handles division by zero)

Computes x / y with protection against division by zero.

## Usage

``` r
div_safe(x, y, eps = .Machine$double.eps)
```

## Arguments

- x:

  Numerator (value object or numeric)

- y:

  Denominator (value object or numeric)

- eps:

  Minimum absolute value for denominator (default: .Machine\$double.eps)

## Value

x / sign(y) \* max(abs(y), eps)
