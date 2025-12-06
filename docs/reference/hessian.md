# Compute Hessian matrix via forward-over-reverse automatic differentiation

Computes the Hessian matrix (matrix of second partial derivatives) of a
scalar function with respect to a vector of parameters. Uses
forward-mode AD on top of reverse-mode AD for efficient, accurate
computation.

## Usage

``` r
hessian(loss_fn, params, value_creator = val)
```

## Arguments

- loss_fn:

  A function that takes a list of parameters and returns a scalar
  loss/objective. The function must use femtograd operations.

- params:

  A list of value objects (the parameters)

- value_creator:

  Function to create value objects (default: val). This allows
  customization for different value types.

## Value

A numeric matrix of dimension (n x n) where n = length(params),
containing the Hessian d²f/dθᵢdθⱼ

## Details

The Hessian is computed using forward-over-reverse mode AD:

1.  For each parameter i, create dual numbers where:

    - Primal = value object holding parameter value

    - Tangent = value object (1 for param i, 0 for others)

2.  Run the loss function with these dual-value inputs

3.  The tangent of the loss is df/dθᵢ (as a value expression)

4.  Call backward() on the tangent loss

5.  Each tangent parameter's gradient gives d²f/dθⱼdθᵢ

This is the "forward-over-reverse" pattern: forward-mode (dual numbers)
computes df/dθᵢ symbolically, then reverse-mode differentiates that to
get d²f/dθⱼdθᵢ.

## Examples

``` r
if (FALSE) { # \dontrun{
# Hessian of f(x,y) = x^2 + x*y + y^2
loss_fn <- function(p) {
  x <- p[[1]]
  y <- p[[2]]
  x^2 + x*y + y^2
}
params <- list(val(1), val(2))
H <- hessian(loss_fn, params)
# H should be [[2, 1], [1, 2]]
} # }
```
