# Constructor for femtofit objects

Creates a fitted model object from MLE results. This is primarily an
internal constructor; users typically obtain femtofit objects from the
[`fit()`](https://queelius.github.io/femtograd/reference/fit.md)
function.

## Usage

``` r
femtofit(
  coefficients,
  vcov,
  loglik,
  hessian = NULL,
  gradient = NULL,
  converged = TRUE,
  iterations = NA_integer_,
  nobs = NULL,
  call = NULL
)

# S3 method for class 'femtofit'
coef(object, ...)

# S3 method for class 'femtofit'
vcov(object, ...)

# S3 method for class 'femtofit'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'femtofit'
logLik(object, ...)

# S3 method for class 'femtofit'
nobs(object, ...)

# S3 method for class 'femtofit'
print(x, ...)

# S3 method for class 'femtofit'
summary(object, ...)

# S3 method for class 'summary.femtofit'
print(x, ...)
```

## Arguments

- coefficients:

  Named numeric vector of parameter estimates.

- vcov:

  Variance-covariance matrix (named).

- loglik:

  Log-likelihood value at the MLE.

- hessian:

  Hessian matrix at the MLE.

- gradient:

  Gradient vector at the MLE (should be near zero).

- converged:

  Logical indicating convergence.

- iterations:

  Number of iterations performed.

- nobs:

  Number of observations (optional).

- call:

  The original function call (optional).

- object:

  A `femtofit` object.

- ...:

  Additional arguments (ignored).

- parm:

  Parameter names or indices (default: all parameters).

- level:

  Confidence level (default: 0.95).

- x:

  A `femtofit` object.

## Value

A `femtofit` object.

## Methods (by generic)

- `coef(femtofit)`: Extract coefficient estimates

- `vcov(femtofit)`: Extract variance-covariance matrix

- `confint(femtofit)`: Compute confidence intervals

- `logLik(femtofit)`: Extract log-likelihood (enables AIC/BIC)

- `nobs(femtofit)`: Number of observations

- `print(femtofit)`: Print fitted model

- `summary(femtofit)`: Summary of fitted model

## Functions

- `print(summary.femtofit)`: Print summary
