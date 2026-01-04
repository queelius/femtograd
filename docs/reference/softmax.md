# Softmax function (numerically stable)

Computes softmax(x)\_i = exp(x_i) / sum(exp(x)) in a numerically stable
way.

## Usage

``` r
softmax(x)
```

## Arguments

- x:

  A value object or numeric vector

## Value

A value object (or numeric) with softmax probabilities

## Details

Uses the identity softmax(x) = softmax(x - max(x)) to prevent overflow.
