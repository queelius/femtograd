# Normal (Gaussian) log-likelihood

Computes the log-likelihood for i.i.d. normal observations. L(μ,σ²\|x) =
-n/2 log(2π) - n/2 log(σ²) - Σ(xᵢ-μ)²/(2σ²)

## Usage

``` r
loglik_normal(mu, sigma, x)
```

## Arguments

- mu:

  Mean parameter (value object)

- sigma:

  Standard deviation parameter (value object), must be positive

- x:

  Numeric vector of observations

## Value

A value object representing the log-likelihood

## Examples

``` r
if (FALSE) { # \dontrun{
x <- rnorm(100, mean = 5, sd = 2)
mu <- val(0)
sigma <- val(1)
ll <- loglik_normal(mu, sigma, x)
backward(ll)
grad(mu)  # score for mu
} # }
```
