# Check if standard errors are reliable

A quick check of whether standard errors from the variance-covariance
matrix are likely to be reliable.

## Usage

``` r
se_reliable(object)
```

## Arguments

- object:

  A `femtofit` object.

## Value

Logical indicating whether SEs are likely reliable.

## Details

Standard errors may be unreliable if:

- The Hessian is not negative definite

- The model did not converge

- Any SE is NA or negative
