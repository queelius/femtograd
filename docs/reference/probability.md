# Transform to probability values

Maps unconstrained values to (0, 1) via the sigmoid function. Use this
for probability parameters.

## Usage

``` r
probability(x)
```

## Arguments

- x:

  Unconstrained value (scalar, vector, or value object)

## Value

Value in (0, 1)

## Details

The transformation is: probability(x) = 1 / (1 + exp(-x)) = sigmoid(x)

For optimization:

      result <- fit(
        function(logit_p) {
          p <- probability(logit_p)
          loglik_bernoulli(p, data)
        },
        params = c(logit_p = 0)  # sigmoid(0) = 0.5
      )
      # To recover p: sigmoid(coef(result)["logit_p"])

## Examples

``` r
if (FALSE) { # \dontrun{
# Parameter that must be in (0, 1)
logit_p <- val(2)
p <- probability(logit_p)  # sigmoid(2) ≈ 0.88
data(p)
} # }
```
