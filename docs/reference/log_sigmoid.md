# Log-sigmoid (numerically stable)

Computes log(sigmoid(x)) = -log(1 + exp(-x)) stably.

## Usage

``` r
log_sigmoid(x)
```

## Arguments

- x:

  A value object or numeric

## Value

log(sigmoid(x))

## Details

Direct computation of log(sigmoid(x)) fails for large negative x
(sigmoid underflows to 0). This uses:

- For x \>= 0: -log(1 + exp(-x)) = -softplus(-x)

- For x \< 0: x - log(1 + exp(x)) = x - softplus(x)
