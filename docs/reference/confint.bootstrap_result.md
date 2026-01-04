# Confidence intervals from bootstrap

Computes confidence intervals using bootstrap methods.

## Usage

``` r
# S3 method for class 'bootstrap_result'
confint(
  object,
  parm = NULL,
  level = 0.95,
  type = c("percentile", "basic", "normal"),
  ...
)
```

## Arguments

- object:

  A `bootstrap_result` object.

- parm:

  Parameter names or indices (default: all).

- level:

  Confidence level (default 0.95).

- type:

  Type of interval: "percentile" (default), "basic", or "normal".

- ...:

  Additional arguments (ignored).

## Value

Matrix of confidence intervals.

## Details

Available methods:

- percentile:

  Uses quantiles of bootstrap distribution directly

- basic:

  Uses 2\*theta_hat - quantiles (pivot method)

- normal:

  Uses bootstrap SE with normal quantiles
