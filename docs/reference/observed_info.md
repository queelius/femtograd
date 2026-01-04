# Observed Fisher information matrix

Returns the observed Fisher information matrix, which is the negative
Hessian of the log-likelihood at the MLE.

## Usage

``` r
observed_info(object, ...)

# S3 method for class 'femtofit'
observed_info(object, ...)
```

## Arguments

- object:

  A fitted model object.

- ...:

  Additional arguments.

## Value

The observed information matrix.

## Methods (by class)

- `observed_info(femtofit)`: Observed information for femtofit
