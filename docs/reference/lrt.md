# Likelihood Ratio Test

Computes the likelihood ratio test (LRT) for comparing nested models.
Works with femtofit objects or raw log-likelihood values.

## Usage

``` r
lrt(null_model, alt_model, df = NULL)
```

## Arguments

- null_model:

  Either a `femtofit` object (simpler model) or a numeric log-likelihood
  value.

- alt_model:

  Either a `femtofit` object (more complex model) or a numeric
  log-likelihood value.

- df:

  Degrees of freedom. If NULL and both arguments are femtofit objects,
  computed as the difference in number of parameters.

## Value

A `likelihood_ratio_test` object (also inherits from `hypothesis_test`)
containing:

- stat:

  The LRT statistic: -2 \* (loglik_null - loglik_alt)

- p.value:

  P-value from chi-squared distribution

- dof:

  Degrees of freedom

- null_loglik:

  Log-likelihood of null model

- alt_loglik:

  Log-likelihood of alternative model

## Details

The likelihood ratio test statistic is: \$\$\Lambda = -2 (\ell_0 -
\ell_1)\$\$ where \\\ell_0\\ is the log-likelihood under the null
(simpler) model and \\\ell_1\\ is the log-likelihood under the
alternative (complex) model.

Under the null hypothesis and regularity conditions, \\\Lambda\\
asymptotically follows a chi-squared distribution with degrees of
freedom equal to the difference in number of free parameters.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate data
set.seed(42)
x <- rnorm(100, mean = 5, sd = 2)

# Fit full model (estimate both mu and sigma)
full <- fit(
  function(mu, log_sigma) {
    sigma <- exp(log_sigma)
    loglik_normal(mu, sigma, x)
  },
  params = c(mu = 0, log_sigma = 0)
)

# Fit null model (mu fixed at 0)
null <- fit(
  function(log_sigma) {
    sigma <- exp(log_sigma)
    loglik_normal(0, sigma, x)  # mu = 0
  },
  params = c(log_sigma = 0)
)

# Likelihood ratio test
test <- lrt(null, full)
test
pval(test)
is_significant_at(test, 0.05)

# Can also use raw log-likelihoods
lrt(null_loglik = -150, alt_loglik = -140, df = 1)
} # }
```
