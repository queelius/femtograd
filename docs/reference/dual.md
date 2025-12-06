# dual R6 class for forward-mode automatic differentiation

dual R6 class for forward-mode automatic differentiation

dual R6 class for forward-mode automatic differentiation

## Details

Represents a dual number (primal, tangent) for computing directional
derivatives via forward-mode AD. When combined with reverse-mode
(value), enables Hessian computation via forward-over-reverse.

Forward-mode AD propagates derivatives alongside values during the
forward pass. For a function f(x), setting tangent(x) = 1 gives f'(x) in
the tangent of the output.

For Hessian computation (forward-over-reverse):

- primal is a `value` object (for reverse-mode gradient)

- tangent is also a `value` object (tracks d/dx of the gradient)

- After backward() on primal, tangent of grad gives Hessian entries

## Public fields

- `primal`:

  The function value (can be numeric or value object)

- `tangent`:

  The directional derivative (can be numeric or value object)

## Methods

### Public methods

- [`dual$new()`](#method-dual-new)

- [`dual$clone()`](#method-dual-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new dual number

#### Usage

    dual$new(primal, tangent = 0)

#### Arguments

- `primal`:

  The primal value (numeric or value object)

- `tangent`:

  The tangent (derivative direction), default 0

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    dual$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
