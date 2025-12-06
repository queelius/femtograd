# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`femtograd` is an R package implementing automatic differentiation (AD) via backpropagation using computational graphs. It's inspired by Karpathy's micrograd but implemented in R using R6 classes. The focus is on **pedagogy and modern statistics** rather than large-scale ML.

Key capabilities:
- First-order gradients via reverse-mode AD (backpropagation)
- Second-order derivatives (Hessians) via forward-over-reverse AD
- Built-in log-likelihood functions for exponential family distributions
- MLE optimization with Newton-Raphson and gradient ascent
- Statistical inference utilities (standard errors, confidence intervals)

## Core Architecture

### Computational Graph System (R/value.R, R/ops.R)

The package uses a directed acyclic graph (DAG) to represent mathematical expressions:

- **value**: R6 class representing nodes. Each node stores:
  - `data`: The numeric value
  - `grad`: Gradient computed during backpropagation
  - `backward_fn`: Function defining local gradient computation
  - `prev`: List of parent nodes

- **Operations**: S3 methods for differentiable operations (+, -, *, /, ^, log, exp, sqrt, sin, cos, lgamma, digamma, sigmoid, etc.)

### Forward-Mode AD for Hessians (R/dual.R, R/hessian.R)

Second derivatives are computed via forward-over-reverse mode:

- **dual**: Class for forward-mode AD. Stores (primal, tangent) pairs.
- **hessian()**: Computes full Hessian matrix by:
  1. Creating dual numbers with value-object primals and tangents
  2. Running forward pass (builds both primal and tangent graphs)
  3. The tangent expression represents df/dθᵢ
  4. Running backward on tangent gives d²f/dθⱼdθᵢ

### Statistical Functions (R/distributions.R, R/optimize.R)

- **Log-likelihoods**: `loglik_normal`, `loglik_exponential`, `loglik_poisson`, `loglik_binomial`, `loglik_gamma`, `loglik_beta`
- **Optimization**: `gradient_ascent`, `gradient_descent`, `newton_raphson`, `find_mle`
- **Inference**: `fisher_information`, `std_errors`, `vcov_matrix`, `confint_mle`, `wald_test`

## Development Commands

### Build and Documentation
```r
devtools::document()  # Generate roxygen2 docs
devtools::build()
devtools::install()
```

### Testing
```r
devtools::test()
covr::package_coverage()
```

Test files:
- `test-ops.R`: Arithmetic, transcendentals, activation functions
- `test-numerical-gradients.R`: Gradient validation via finite differences
- `test-math-ops.R`: lgamma, digamma, sin, cos, etc.
- `test-dual.R`: Forward-mode AD and dual number operations
- `test-hessian.R`: Hessian computation
- `test-distributions.R`: Log-likelihood functions
- `test-optimize.R`: Optimization and MLE routines

## Key Design Principles

1. **Pedagogy over performance**: Code is written for clarity, not speed
2. **Statistics focus**: Designed for MLE, Hessians, and inference rather than deep learning
3. **Graph Construction**: Operations automatically build computational graphs
4. **Lazy Gradient Computation**: Gradients only computed on `backward()` call
5. **Gradient Accumulation**: Parent node gradients accumulate from multiple children

## Important Implementation Details

### Operator Dispatch Limitation
R dispatches S3 methods based on the **first** argument. This means:
```r
x <- val(5)
x * 2     # Works: dispatches to *.value
2 * x     # FAILS: dispatches to *.default (numeric)
```
Always put the `value` object first in binary operations.

### Gradient Update Pattern
```r
# WRONG: Creates new value, breaks graph
rate <- rate + lr * grad(rate)

# CORRECT: Update data in place
data(rate) <- data(rate) + lr * grad(rate)
```

### Hessian Computation
```r
loss_fn <- function(p) {
  mu <- p[[1]]
  sigma <- p[[2]]
  loglik_normal(mu, sigma, data)
}
params <- list(val(0), val(1))
H <- hessian(loss_fn, params)
se <- std_errors(H, is_loglik = TRUE)
```

### Finding MLE with Standard Errors
```r
result <- find_mle(loglik_fn, initial_params)
result$estimate  # MLE values
result$se        # Standard errors
result$vcov      # Variance-covariance matrix
```
