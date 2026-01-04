# Fit a model via maximum likelihood

Finds the maximum likelihood estimates of parameters and returns a
fitted model object with standard errors, confidence intervals, and
other inference quantities computed automatically via autodiff.

## Usage

``` r
fit(loglik, params, method = c("bfgs", "lbfgs", "newton", "gradient"), ...)
```

## Arguments

- loglik:

  A log-likelihood function. Can be specified in two ways:

  - Named arguments: `function(mu, sigma) loglik_normal(mu, sigma, x)`

  - Single parameter argument:
    `function(p) loglik_normal(p$mu, p$sigma, x)`

  The function should return a scalar value object.

- params:

  Named numeric vector of initial parameter values, e.g.,
  `c(mu = 0, sigma = 1)`.

- method:

  Optimization method: "bfgs" (default), "lbfgs", "newton", or
  "gradient".

- ...:

  Additional arguments passed to the optimizer.

## Value

A `femtofit` object containing:

- Coefficient estimates (accessible via
  [`coef()`](https://rdrr.io/r/stats/coef.html))

- Variance-covariance matrix (accessible via
  [`vcov()`](https://rdrr.io/r/stats/vcov.html))

- Log-likelihood value (accessible via
  [`logLik()`](https://rdrr.io/r/stats/logLik.html))

- Standard errors, Hessian, convergence info

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate data
set.seed(42)
x <- rnorm(100, mean = 5, sd = 2)

# Fit using named arguments (recommended)
result <- fit(
  function(mu, sigma) loglik_normal(mu, sigma, x),
  params = c(mu = 0, sigma = 1)
)

# Or using single parameter argument
result <- fit(
  function(p) loglik_normal(p$mu, p$sigma, x),
  params = c(mu = 0, sigma = 1)
)

# Standard R generics work
coef(result)      # Parameter estimates
vcov(result)      # Variance-covariance matrix
confint(result)   # 95% confidence intervals
logLik(result)    # Log-likelihood (works with AIC/BIC)
AIC(result)       # Akaike information criterion
summary(result)   # Full summary with p-values
} # }
```
