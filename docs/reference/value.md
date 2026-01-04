# value R6 class

value R6 class

value R6 class

## Details

Represents a node in the computational graph with automatic
differentiation. Each value object holds data as a matrix, gradient as a
matrix of the same dimensions, a backward function, previous nodes, and
an optional label.

All data in femtograd uses matrix representation:

- Scalars are 1x1 matrices

- Column vectors are n x 1 matrices

- Row vectors are 1 x n matrices

- Matrices are m x n matrices

This ensures consistent behavior for data(), grad(), and hessian().

## Public fields

- `data`:

  Numeric matrix containing the value

- `grad`:

  Gradient matrix (same dimensions as data), initially zeros

- `backward_fn`:

  A function that performs the backward pass (gradient computation)

- `prev`:

  A list of previous nodes in the computational graph

- `label`:

  Optional character label for debugging (default: "")

## Methods

### Public methods

- [`value$new()`](#method-value-new)

- [`value$clone()`](#method-value-clone)

------------------------------------------------------------------------

### Method `new()`

Initializes a new value object with the given data, list of children,
and optional label.

#### Usage

    value$new(data, children = list(), label = "")

#### Arguments

- `data`:

  Numeric value (scalar, vector, or matrix) - will be converted to
  matrix

- `children`:

  List of previous nodes in the computational graph (default: empty
  list)

- `label`:

  Optional character label for debugging

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    value$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
