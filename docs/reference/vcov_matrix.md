# Compute variance-covariance matrix from Hessian

Compute variance-covariance matrix from Hessian

## Usage

``` r
vcov_matrix(hess, is_loglik = TRUE)
```

## Arguments

- hess:

  The Hessian matrix

- is_loglik:

  If TRUE, returns -H⁻¹ (for log-likelihood)

## Value

The variance-covariance matrix
