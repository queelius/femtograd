# value R6 class

value R6 class

value R6 class

## Details

Represents a node in the computational graph with automatic
differentiation. Each value object holds data, gradient, a backward
function, previous nodes, and an optional label for debugging.

## Public fields

- `data`:

  Numeric value of the object

- `grad`:

  Gradient value, initially set to 0

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

  Numeric value of the object

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
