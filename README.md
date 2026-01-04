# femtograd

Automatic differentiation for statistical computing in R.

`femtograd` provides reverse-mode AD (backpropagation) for first-order gradients and forward-over-reverse AD for Hessian computation. Designed for **pedagogy and modern statistics** rather than large-scale ML.

## Features

- **Reverse-mode AD**: Efficient gradient computation via backpropagation
- **Forward-over-reverse AD**: Second-order derivatives (Hessians) for Newton-Raphson and Fisher information
- **Log-likelihood functions**: Exponential family distributions (normal, exponential, Poisson, binomial, gamma, beta, negative binomial, logistic, Bernoulli)
- **MLE optimization**: Gradient ascent/descent, Newton-Raphson, BFGS, L-BFGS with line search
- **Statistical inference**: Standard errors, variance-covariance matrices, confidence intervals, Wald tests
- **Numerical stability**: logsumexp, softmax, safe log/division, overflow protection

## Installation

```r
# Install from GitHub
devtools::install_github("queelius/femtograd")
```

## Quick Start

```r
library(femtograd)

# Generate exponential data
set.seed(42)
x <- rexp(100, rate = 2)

# Define log-likelihood
loglik <- function(p) loglik_exponential(p[[1]], x)

# Find MLE with standard errors
result <- find_mle(loglik, list(val(1)), method = "newton")
result$estimate  # MLE
result$se        # Standard errors
```

## Core Concepts

### Computational Graph

Create differentiable values with `val()`:

```r
x <- val(3)
y <- val(4)
z <- x^2 + x*y  # z = 9 + 12 = 21

backward(z)
grad(x)  # dz/dx = 2x + y = 10
grad(y)  # dz/dy = x = 3
```

### Supported Operations

Arithmetic: `+`, `-`, `*`, `/`, `^`

Transcendental: `log`, `exp`, `sqrt`, `sin`, `cos`, `tanh`, `lgamma`, `digamma`, `log1p`

Activations: `sigmoid`, `relu`, `softplus`, `logit`

Stability: `logsumexp`, `softmax`, `log_safe`, `div_safe`, `exp_safe`, `sigmoid_stable`, `log_sigmoid`

### Hessian Computation

```r
f <- function(p) {
  x <- p[[1]]
  y <- p[[2]]
  x^2 + x*y + y^2
}

H <- hessian(f, list(val(1), val(2)))
# Returns 2x2 Hessian matrix
```

## Maximum Likelihood Estimation

### Built-in Distributions

```r
# Normal distribution
loglik <- function(p) loglik_normal(p[[1]], p[[2]], data)

# Poisson distribution
loglik <- function(p) loglik_poisson(p[[1]], data)

# Gamma distribution
loglik <- function(p) loglik_gamma(p[[1]], p[[2]], data)
```

### Optimization Methods

```r
# Gradient ascent
result <- gradient_ascent(loglik, params, lr = 0.01, max_iter = 1000)

# Newton-Raphson (uses Hessian)
result <- newton_raphson(loglik, params, max_iter = 50)

# BFGS quasi-Newton
result <- bfgs(loglik, params, maximize = TRUE)

# L-BFGS (memory efficient)
result <- lbfgs(loglik, params, m = 10, maximize = TRUE)

# Full MLE with inference
result <- find_mle(loglik, params, method = "newton")
```

### Statistical Inference

```r
result <- find_mle(loglik, params)

# Standard errors from Fisher information
result$se

# Variance-covariance matrix
result$vcov

# Confidence intervals
confint_mle(result, level = 0.95)

# Wald test
wald_test(result, null_values = c(0, 1))
```

## Example: Normal MLE

```r
set.seed(123)
true_mu <- 5
true_sigma <- 2
x <- rnorm(100, true_mu, true_sigma)

loglik <- function(p) loglik_normal(p[[1]], p[[2]], x)

# Find MLE
result <- find_mle(loglik, list(val(mean(x)), val(sd(x))))

cat("MLE mu:", result$estimate[1], "SE:", result$se[1], "\n")
cat("MLE sigma:", result$estimate[2], "SE:", result$se[2], "\n")

# 95% confidence intervals
confint_mle(result)
```

## Design Philosophy

- **Pedagogy over performance**: Code prioritizes clarity for teaching AD concepts
- **Statistics focus**: Built for MLE, Hessians, and inference rather than deep learning
- **Composable foundation**: Designed as a building block for other statistical packages

## Acknowledgments

Inspired by Karpathy's [micrograd](https://github.com/karpathy/micrograd).

## License

GPL-3
