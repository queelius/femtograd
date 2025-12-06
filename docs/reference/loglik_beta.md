# Beta distribution log-likelihood

Computes the log-likelihood for i.i.d. beta observations. L(α,β\|x) =
n\*\[log(Γ(α+β)) - log(Γ(α)) - log(Γ(β))\] + (α-1)\*Σlog(xᵢ) +
(β-1)\*Σlog(1-xᵢ)

## Usage

``` r
loglik_beta(alpha, beta, x)
```

## Arguments

- alpha:

  Shape parameter α (value object), must be positive

- beta:

  Shape parameter β (value object), must be positive

- x:

  Numeric vector of observations in (0,1)

## Value

A value object representing the log-likelihood
