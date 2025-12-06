# Log(1+x) for value objects

Computes log(1+x) in a numerically stable way. Essential for computing
log-likelihoods with probabilities near 0 or 1.

## Usage

``` r
# S3 method for class 'value'
log1p(x)
```

## Arguments

- x:

  A value object

## Value

A new value object representing log(1+x)
