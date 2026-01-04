# Check convergence diagnostics

Examines convergence properties of a fitted model.

## Usage

``` r
check_convergence(object, tol = 1e-04)
```

## Arguments

- object:

  A `femtofit` object.

- tol:

  Tolerance for gradient norm (default 1e-4).

## Value

A `convergence_check` object containing:

- converged:

  From the optimizer

- gradient_norm:

  Norm of gradient at solution

- gradient_ok:

  Is gradient norm below tolerance?

- iterations:

  Number of iterations used

- loglik:

  Final log-likelihood value

- warnings:

  Character vector of any issues detected

## Details

At a maximum, the gradient should be (near) zero. Large gradient norms
may indicate premature termination or numerical issues.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
x <- rnorm(100, mean = 5, sd = 2)

result <- fit(
  function(mu, log_sigma) loglik_normal(mu, exp(log_sigma), x),
  params = c(mu = 0, log_sigma = 0)
)

check <- check_convergence(result)
check
} # }
```
