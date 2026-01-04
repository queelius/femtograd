# Pareto distribution log-likelihood

Computes the log-likelihood for i.i.d. Pareto observations. L(α, xₘ \|
x) = n*log(α) + n*α\*log(xₘ) - (α+1)\*Σlog(xᵢ)

## Usage

``` r
loglik_pareto(alpha, x_min, x)
```

## Arguments

- alpha:

  Shape parameter α (value object), must be positive

- x_min:

  Minimum/scale parameter xₘ (fixed positive number). All observations
  must be \>= x_min.

- x:

  Numeric vector of observations (must be \>= x_min)

## Value

A value object representing the log-likelihood

## Details

The Pareto distribution is used to model heavy-tailed phenomena like
income distributions, city sizes, etc. Here x_min is typically known
(e.g., min(x)) and alpha is estimated.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate Pareto data
alpha_true <- 2
x_min <- 1
u <- runif(100)
x <- x_min * (1 - u)^(-1/alpha_true)

# Fit (alpha only, x_min = min(x) is fixed)
result <- fit(
  function(log_alpha) {
    alpha <- exp(log_alpha)
    loglik_pareto(alpha, x_min = min(x), x)
  },
  params = c(log_alpha = 0)
)
} # }
```
