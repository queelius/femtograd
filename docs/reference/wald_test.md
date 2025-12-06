# Wald test for hypothesis testing

Computes Wald test statistic for H0: θ = θ₀ W = (θ̂ - θ₀)' I(θ̂) (θ̂ - θ₀)
~ χ²(k)

## Usage

``` r
wald_test(mle_result, null_values = NULL)
```

## Arguments

- mle_result:

  Result from find_mle

- null_values:

  Null hypothesis values (default: 0 for all params)

## Value

A list with test statistic, df, and p-value
