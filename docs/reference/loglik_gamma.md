# Gamma distribution log-likelihood

Computes the log-likelihood for i.i.d. gamma observations. L(α,β\|x) =
n*α*log(β) - n\*log(Γ(α)) + (α-1)*Σlog(xᵢ) - β*Σxᵢ

## Usage

``` r
loglik_gamma(shape, rate, x)
```

## Arguments

- shape:

  Shape parameter α (value object), must be positive

- rate:

  Rate parameter β (value object), must be positive

- x:

  Numeric vector of observations (must be positive)

## Value

A value object representing the log-likelihood

## Details

The gamma distribution is parameterized with shape α and rate β, where
EX = α/β and VarX = α/β².
