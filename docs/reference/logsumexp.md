# Log-Sum-Exp (numerically stable)

Computes log(sum(exp(x))) in a numerically stable way by factoring out
the maximum value: log(sum(exp(x))) = m + log(sum(exp(x - m))) where m =
max(x).

## Usage

``` r
logsumexp(..., na.rm = FALSE)
```

## Arguments

- ...:

  value objects and/or numeric vectors to include in the sum

- na.rm:

  Logical, whether to remove NA values

## Value

A value object representing log(sum(exp(x)))

## Details

The naive computation of log(sum(exp(x))) overflows when x contains
large values (exp(710) = Inf in double precision) and underflows when x
contains very negative values. This implementation is stable for any
finite input.

## Examples

``` r
if (FALSE) { # \dontrun{
x <- val(c(1000, 1001, 1002))  # Would overflow with naive exp()
logsumexp(x)  # Returns ~1002.41 correctly
} # }
```
