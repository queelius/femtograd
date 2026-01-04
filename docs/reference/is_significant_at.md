# Check if test is significant at given level

Check if test is significant at given level

## Usage

``` r
is_significant_at(x, alpha = 0.05, ...)

# S3 method for class 'hypothesis_test'
is_significant_at(x, alpha = 0.05, ...)
```

## Arguments

- x:

  A hypothesis test object

- alpha:

  Significance level (default 0.05)

- ...:

  Additional arguments (ignored)

## Value

Logical indicating significance

## Methods (by class)

- `is_significant_at(hypothesis_test)`: Significance check for
  hypothesis tests
