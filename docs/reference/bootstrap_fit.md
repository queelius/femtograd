# Bootstrap standard errors and confidence intervals

Performs bootstrap resampling to estimate the sampling distribution of
parameter estimates.

## Usage

``` r
bootstrap_fit(loglik_maker, data, params, n_boot = 500, progress = TRUE, ...)
```

## Arguments

- loglik_maker:

  A function that takes data and returns a log-likelihood function
  suitable for
  [`fit()`](https://queelius.github.io/femtograd/reference/fit.md). See
  examples.

- data:

  The data (vector, matrix, or data frame).

- params:

  Initial parameter values for
  [`fit()`](https://queelius.github.io/femtograd/reference/fit.md).

- n_boot:

  Number of bootstrap replications (default 500).

- progress:

  Logical; show progress messages (default TRUE).

- ...:

  Additional arguments passed to
  [`fit()`](https://queelius.github.io/femtograd/reference/fit.md).

## Value

A `bootstrap_result` object containing:

- estimates:

  Matrix of bootstrap estimates (n_boot x n_params)

- original:

  Original parameter estimates

- se:

  Bootstrap standard errors

- bias:

  Estimated bias

## Details

The `loglik_maker` function should accept data and return a
log-likelihood function that can be passed to
[`fit()`](https://queelius.github.io/femtograd/reference/fit.md). This
design allows bootstrap to work with any model specification.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
x <- rnorm(50, mean = 5, sd = 2)

# Define a function that creates the loglik function from data
make_loglik <- function(data) {
  function(mu, log_sigma) {
    loglik_normal(mu, exp(log_sigma), data)
  }
}

# Bootstrap
boot <- bootstrap_fit(
  make_loglik,
  data = x,
  params = c(mu = 0, log_sigma = 0),
  n_boot = 200
)

# Results
boot$se          # Bootstrap SEs
confint(boot)    # Bootstrap CIs

# Compare with Wald SEs
result <- fit(make_loglik(x), params = c(mu = 0, log_sigma = 0))
se(result)       # Wald SEs
} # }
```
