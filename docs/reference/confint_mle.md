# Compute confidence intervals from MLE results

Compute confidence intervals from MLE results

## Usage

``` r
confint_mle(mle_result, level = 0.95, param_names = NULL)
```

## Arguments

- mle_result:

  Result from find_mle

- level:

  Confidence level, default 0.95

- param_names:

  Optional names for parameters

## Value

A matrix with columns: estimate, se, lower, upper
