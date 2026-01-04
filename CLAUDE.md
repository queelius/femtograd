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
  - `data`: Numeric matrix (all values are matrices internally)
  - `grad`: Gradient matrix (same dimensions as data)
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

### Statistical Functions (R/distributions.R, R/optimize.R, R/fit.R)

- **Log-likelihoods**: `loglik_normal`, `loglik_exponential`, `loglik_poisson`, `loglik_binomial`, `loglik_gamma`, `loglik_beta`
- **Optimization**: `gradient_ascent`, `gradient_descent`, `newton_raphson`, `bfgs`, `lbfgs`
- **Model fitting**: `fit()` returns `femtofit` object with standard R generics
- **Inference**: `coef()`, `vcov()`, `confint()`, `logLik()`, `AIC()`, `BIC()`, `summary()`

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
- `test-fit.R`: fit() function and femtofit class with R generics
- `test-analytical.R`: Known gradients/Hessians (power rule, chain rule, etc.)
- `test-stability.R`: Numerical stability utilities (logsumexp, softmax, etc.)

## Key Design Principles

1. **Pedagogy over performance**: Code is written for clarity, not speed
2. **Statistics focus**: Designed for MLE, Hessians, and inference rather than deep learning
3. **Base R generics**: Uses standard R generics (`coef`, `vcov`, `confint`, `logLik`) for interoperability
4. **Matrix representation**: All values are matrices internally (see DESIGN.md)
5. **Graph Construction**: Operations automatically build computational graphs
6. **Lazy Gradient Computation**: Gradients only computed on `backward()` call
7. **Gradient Accumulation**: Parent node gradients accumulate from multiple children

## Important Implementation Details

### Matrix Representation
All data is stored as matrices (see DESIGN.md for rationale):
- Scalars are 1×1 matrices
- Vectors are n×1 column matrices
- Matrices are m×n matrices

```r
val(5)           # 1×1 matrix internally
val(c(1,2,3))    # 3×1 column matrix internally
data(x)          # Returns scalar for 1×1 (drop=TRUE default)
data(x, drop=FALSE)  # Always returns matrix
```

### Broadcasting
Scalar (1×1) matrices broadcast over larger matrices:
```r
a <- val(2)           # 1×1
x <- val(c(1,2,3))    # 3×1
y <- a * x            # 3×1: c(2,4,6)
backward(y)
grad(a)  # Scalar: sum of gradients = 6
grad(x)  # 3×1: c(2,2,2)
```

### Operator Dispatch
R dispatches S3 methods based on the **first** argument. However, femtograd now handles both orderings:
```r
x <- val(5)
x * 2     # Works: dispatches to *.value
2 * x     # Also works: scalar is converted to value internally
```

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

### Finding MLE with Standard Errors (Recommended)
```r
# Define log-likelihood with named parameters
# Use log-parameterization for constrained parameters (e.g., sigma > 0)
result <- fit(
  function(mu, log_sigma) {
    sigma <- exp(log_sigma)
    loglik_normal(mu, sigma, data)
  },
  params = c(mu = 0, log_sigma = 0)
)

# Standard R generics
coef(result)      # Named parameter estimates
vcov(result)      # Variance-covariance matrix
confint(result)   # 95% confidence intervals
se(result)        # Standard errors
logLik(result)    # Log-likelihood (works with AIC/BIC)
AIC(result)       # Akaike information criterion
summary(result)   # Full summary with z-values and p-values
```

### Alternative: Low-level find_mle()
```r
result <- find_mle(loglik_fn, initial_params)
result$estimate  # MLE values
result$se        # Standard errors
result$vcov      # Variance-covariance matrix
```

### Parameter Constraints
For constrained parameters, use transformations:
```r
# sigma > 0: use log_sigma, then exp(log_sigma) in likelihood
# 0 < p < 1: use logit_p, then sigmoid(logit_p) in likelihood

result <- fit(
  function(mu, log_sigma) {
    sigma <- exp(log_sigma)  # Ensures sigma > 0
    loglik_normal(mu, sigma, data)
  },
  params = c(mu = 0, log_sigma = 0)
)

# To get sigma estimate: exp(coef(result)["log_sigma"])
```
