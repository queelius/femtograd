# Logistic regression log-likelihood (binary)

Computes the log-likelihood for binary logistic regression. L(β\|X,y) =
Σyᵢ\*log(pᵢ) + (1-yᵢ)\*log(1-pᵢ) where pᵢ = sigmoid(Xᵢ·β)

## Usage

``` r
loglik_logistic(beta, X, y)
```

## Arguments

- beta:

  Coefficient vector (list of value objects)

- X:

  Design matrix (n x p numeric matrix)

- y:

  Binary response vector (0 or 1)

## Value

A value object representing the log-likelihood

## Details

Uses the numerically stable form: log(p) = -log(1 + exp(-η)) and
log(1-p) = -log(1 + exp(η)) where η = Xβ
