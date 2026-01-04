# Backtracking line search (Armijo condition)

Finds a step size that satisfies the Armijo condition for sufficient
decrease.

## Usage

``` r
line_search(
  objective_fn,
  param_values,
  direction,
  grad,
  maximize = FALSE,
  alpha = 1,
  c1 = 1e-04,
  rho = 0.5,
  max_iter = 20
)
```

## Arguments

- objective_fn:

  Function to minimize/maximize

- param_values:

  Current parameter values (numeric vector)

- direction:

  Search direction (numeric vector)

- grad:

  Current gradient (numeric vector)

- maximize:

  If TRUE, maximize; if FALSE, minimize

- alpha:

  Initial step size, default 1

- c1:

  Armijo parameter (sufficient decrease), default 1e-4

- rho:

  Backtracking factor, default 0.5

- max_iter:

  Maximum backtracking iterations, default 20

## Value

Step size satisfying Armijo condition

## Details

Armijo condition (for minimization): f(x + α*d) \<= f(x) + c1*α\*g'\*d
where d is the search direction and g is the gradient.
