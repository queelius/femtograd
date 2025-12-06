# Compute standard errors from Hessian

Extracts standard errors from the Hessian of a log-likelihood. SE(θ̂ᵢ) =
sqrt(I(θ̂)⁻¹ᵢᵢ) = sqrt(-H⁻¹ᵢᵢ)

## Usage

``` r
std_errors(hess, is_loglik = TRUE)
```

## Arguments

- hess:

  The Hessian matrix (from
  [`hessian()`](https://queelius.github.io/femtograd/reference/hessian.md))

- is_loglik:

  If TRUE, treats hess as Hessian of log-likelihood and uses -H⁻¹. If
  FALSE, uses H⁻¹ directly.

## Value

A numeric vector of standard errors
