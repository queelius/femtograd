# Wald Test for Model Parameters

Performs Wald tests for individual parameters or linear combinations of
parameters in a fitted model.

## Usage

``` r
wald_test(object, ...)

# S3 method for class 'femtofit'
wald_test(object, parm = NULL, null_value = 0, ...)

# Default S3 method
wald_test(object, se, null_value = 0, ...)
```

## Arguments

- object:

  A `femtofit` object, or for the generic: the parameter estimate.

- ...:

  Additional arguments.

- parm:

  Parameter name(s) or indices to test. If NULL, tests all parameters.

- null_value:

  The null hypothesis value(s). Default is 0.

- se:

  Standard error (for default method)

## Value

For single parameter tests, a `wald_test` object (also inherits from
`hypothesis_test`) containing:

- stat:

  The Wald chi-squared statistic

- p.value:

  P-value from chi-squared(1) distribution

- dof:

  Degrees of freedom (1 for single parameter)

- z:

  The z-score

- estimate:

  Parameter estimate

- se:

  Standard error

- null_value:

  The null hypothesis value

For multiple parameters, returns a list of wald_test objects.

## Details

The Wald test statistic for a single parameter is: \$\$W =
\left(\frac{\hat{\theta} - \theta_0}{SE(\hat{\theta})}\right)^2\$\$

which follows a chi-squared distribution with 1 degree of freedom under
H0.

## Methods (by class)

- `wald_test(femtofit)`: Wald test for femtofit objects

- `wald_test(default)`: Wald test from raw estimate and SE

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
x <- rnorm(100, mean = 5, sd = 2)

result <- fit(
  function(mu, log_sigma) loglik_normal(mu, exp(log_sigma), x),
  params = c(mu = 0, log_sigma = 0)
)

# Test if mu = 0
test <- wald_test(result, "mu")
pval(test)

# Test if mu = 5
wald_test(result, "mu", null_value = 5)

# Test all parameters
wald_test(result)
} # }
```
