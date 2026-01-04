# Check Hessian properties

Examines the Hessian matrix for numerical issues that may indicate
problems with the optimization or model specification.

## Usage

``` r
check_hessian(object)
```

## Arguments

- object:

  A `femtofit` object.

## Value

A `hessian_check` object containing:

- is_negative_definite:

  Logical: is -H positive definite? (required for MLE)

- eigenvalues:

  Eigenvalues of -H (observed information)

- condition_number:

  Ratio of largest to smallest eigenvalue

- rank:

  Numerical rank of the matrix

- is_singular:

  Logical: is the matrix numerically singular?

- warnings:

  Character vector of any issues detected

## Details

For a proper maximum of the log-likelihood, the Hessian should be
negative definite (equivalently, -H should be positive definite).

Common issues:

- **Not negative definite**: May indicate a saddle point rather than
  maximum

- **High condition number**: Indicates near-singularity, poorly
  identified parameters

- **Singular matrix**: Parameters are not identifiable from the data

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
x <- rnorm(100, mean = 5, sd = 2)

result <- fit(
  function(mu, log_sigma) loglik_normal(mu, exp(log_sigma), x),
  params = c(mu = 0, log_sigma = 0)
)

check <- check_hessian(result)
check

# Access specific properties
check$is_negative_definite
check$condition_number
} # }
```
