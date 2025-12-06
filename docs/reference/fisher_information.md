# Compute observed Fisher information matrix

For a log-likelihood function, the observed Fisher information is the
negative Hessian evaluated at the MLE. This is used for computing
standard errors of MLEs.

## Usage

``` r
fisher_information(loglik_fn, params)
```

## Arguments

- loglik_fn:

  Log-likelihood function taking a list of value parameters

- params:

  List of value objects at MLE

## Value

The observed Fisher information matrix (negative Hessian)

## Details

The observed Fisher information I(θ̂) = -H(θ̂) where H is the Hessian of
the log-likelihood. Standard errors are sqrt(diag(I⁻¹)).
