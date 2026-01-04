# Compute profile likelihood for a parameter

Computes the profile log-likelihood for a single parameter by maximizing
over all other parameters at each fixed value.

## Usage

``` r
profile_loglik(object, parm, values = NULL, n_points = 20, range_mult = 3, ...)
```

## Arguments

- object:

  A `femtofit` object from
  [`fit()`](https://queelius.github.io/femtograd/reference/fit.md).

- parm:

  The parameter name or index to profile.

- values:

  Optional vector of values at which to evaluate the profile. If NULL, a
  grid is computed based on the MLE and standard error.

- n_points:

  Number of points in the grid (default 20). Ignored if `values` is
  provided.

- range_mult:

  Multiplier for the range around MLE (default 3). The grid spans MLE ±
  range_mult \* SE.

- ...:

  Additional arguments passed to the optimizer.

## Value

A `profile_likelihood` object containing:

- parameter:

  Name of the profiled parameter

- values:

  Grid of parameter values

- profile_loglik:

  Profile log-likelihood at each value

- mle:

  MLE value of the parameter

- max_loglik:

  Maximum log-likelihood

## Details

The profile log-likelihood for parameter θᵢ is defined as:
\$\$pl(\theta_i) = \max\_{\theta\_{-i}} \ell(\theta_i, \theta\_{-i})\$\$

where θ₋ᵢ denotes all parameters except θᵢ.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
x <- rnorm(100, mean = 5, sd = 2)

result <- fit(
  function(mu, log_sigma) loglik_normal(mu, exp(log_sigma), x),
  params = c(mu = 0, log_sigma = 0)
)

# Profile likelihood for mu
prof <- profile_loglik(result, "mu")
plot(prof)
} # }
```
