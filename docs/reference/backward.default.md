# Default implementation does not propagate gradients. For instance, if we have a constant, then the partial of the constant is not meaningful.

Default implementation does not propagate gradients. For instance, if we
have a constant, then the partial of the constant is not meaningful.

## Usage

``` r
# Default S3 method
backward(x, ...)
```

## Arguments

- ...:

  pass additional arguments

- e:

  A value object for which the backward pass should be performed
