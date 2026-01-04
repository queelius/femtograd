# Profile confidence intervals

Computes confidence intervals based on the profile likelihood. These are
more accurate than Wald intervals when the likelihood is non-quadratic.

## Usage

``` r
confint_profile(object, parm = NULL, level = 0.95, ...)
```

## Arguments

- object:

  A `femtofit` object.

- parm:

  Parameter name(s) or indices. If NULL, computes for all.

- level:

  Confidence level (default 0.95).

- ...:

  Additional arguments passed to `profile_loglik`.

## Value

A matrix with columns for lower and upper bounds.

## Details

The profile confidence interval at level 1-α is: \$\$CI = \\\theta_i :
2(\ell(\hat{\theta}) - pl(\theta_i)) \leq \chi^2\_{1,\alpha}\\\$\$

This is the set of parameter values that would not be rejected by a
likelihood ratio test at level α.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
x <- rnorm(100, mean = 5, sd = 2)

result <- fit(
  function(mu, log_sigma) loglik_normal(mu, exp(log_sigma), x),
  params = c(mu = 0, log_sigma = 0)
)

# Profile-based CIs (more accurate for non-quadratic likelihoods)
confint_profile(result)

# Compare with Wald CIs
confint(result)
} # }
```
