#' value R6 class
#'
#' Represents a node in the computational graph with automatic differentiation.
#' Each value object holds data as a matrix, gradient as a matrix of the same
#' dimensions, a backward function, previous nodes, and an optional label.
#'
#' All data in femtograd uses matrix representation:
#' - Scalars are 1x1 matrices
#' - Column vectors are n x 1 matrices
#' - Row vectors are 1 x n matrices
#' - Matrices are m x n matrices
#'
#' This ensures consistent behavior for data(), grad(), and hessian().
#'
#' @field data Numeric matrix containing the value
#' @field grad Gradient matrix (same dimensions as data), initially zeros
#' @field backward_fn A function that performs the backward pass (gradient computation)
#' @field prev A list of previous nodes in the computational graph
#' @field label Optional character label for debugging (default: "")
#'
#' @import R6
#' @export
value <- R6Class("value",
   public = list(
     data = NULL,
     grad = NULL,
     backward_fn = NULL,
     prev = NULL,
     label = "",

     #' @description Initializes a new value object with the given data, list of children, and optional label.
     #' @param data Numeric value (scalar, vector, or matrix) - will be converted to matrix
     #' @param children List of previous nodes in the computational graph (default: empty list)
     #' @param label Optional character label for debugging
     initialize = function(data, children = list(), label = "")
     {
       self$data <- as_matrix(data)
       self$grad <- matrix(0, nrow = nrow(self$data), ncol = ncol(self$data))
       self$backward_fn <- function() NULL
       self$prev <- children
       self$label <- label
     })
)

#' Convert input to matrix representation
#'
#' Internal helper to ensure all data is stored as matrices.
#' - Scalars become 1x1 matrices
#' - Vectors become n x 1 column matrices
#' - Matrices remain unchanged
#'
#' @param x Numeric scalar, vector, or matrix
#' @return A matrix
#' @keywords internal
as_matrix <- function(x)
{
  if (is.matrix(x)) {
    return(x)
  } else if (is.numeric(x)) {
    # Convert scalar or vector to column matrix
    return(matrix(x, ncol = 1))
  } else {
    stop("as_matrix: expected numeric scalar, vector, or matrix")
  }
}

#' Check if matrix is a scalar (1x1)
#' @keywords internal
is_scalar_matrix <- function(x)
{
  is.matrix(x) && nrow(x) == 1 && ncol(x) == 1
}

#' Broadcast matrices for element-wise operations
#'
#' Handles broadcasting when one operand is a 1x1 scalar matrix.
#' Returns the scalar value if broadcasting is needed, otherwise
#' returns the original matrix.
#'
#' @param m1 First matrix
#' @param m2 Second matrix
#' @return List with possibly broadcast matrices
#' @keywords internal
broadcast_matrices <- function(m1, m2)
{
  # If both are same dimensions, no broadcasting needed
  if (identical(dim(m1), dim(m2))) {
    return(list(m1 = m1, m2 = m2))
  }

  # If m1 is 1x1 scalar, extract value for broadcasting

  if (is_scalar_matrix(m1)) {
    return(list(m1 = m1[1,1], m2 = m2))
  }

  # If m2 is 1x1 scalar, extract value for broadcasting
  if (is_scalar_matrix(m2)) {
    return(list(m1 = m1, m2 = m2[1,1]))
  }

  # Otherwise, dimensions must match (R will error on operation)
  list(m1 = m1, m2 = m2)
}

#' \code{value} object constructor
#'
#' Creates a new value node in the computational graph with automatic
#' differentiation. Input is converted to matrix representation internally.
#'
#' @param data Numeric value (scalar, vector, or matrix)
#' @param label Optional character label for debugging
#'
#' @details
#' All values are stored as matrices internally:
#' - `val(5)` creates a 1x1 matrix
#' - `val(c(1,2,3))` creates a 3x1 column vector
#' - `val(matrix(...))` preserves matrix dimensions
#'
#' @return A value object
#' @export
val <- function(data, label = "")
{
  value$new(data, label = label)
}

#' Generic function for the Backward pass for automatic differentiation (finds
#' the gradient of every sub-node in the computational graph with respect to
#' \code{e}). In other words, it is responsible for computing the gradient with
#' respect to \code{e}.
#'
#' This generic function should be implemented for specific classes. We provide
#' an implementation for \code{value} objects.
#'
#' @param e An object for which the backward pass should be performed
#' @param ... additional arguments to pass
#' @export
backward <- function(e, ...)
{
  UseMethod("backward")
}

#' Backward pass for value objects
#'
#' Performs the backward pass (gradient computation) for a value object in the
#' computational graph with automatic differentiation. Initializes the gradient
#' of the output node to a matrix of ones (same shape as data).
#'
#' @param x A value object for which the backward pass should be performed
#' @param ... pass additional arguments
#' @export
backward.value <- function(x, ...)
{
  topo <- topological_sort(x)
  # Initialize gradient to matrix of ones (same shape as data)
  x$grad <- matrix(1, nrow = nrow(x$data), ncol = ncol(x$data))
  for (v in rev(topo))
    v$backward_fn()
}

#' @export
zero_grad <- function(e)
{
  UseMethod("zero_grad")
}

#' @export
zero_grad.value <- function(e)
{
  reset_grad <- function(v) {
    v$grad <- matrix(0, nrow = nrow(v$data), ncol = ncol(v$data))
    for (child in v$prev)
      reset_grad(child)
  }
  reset_grad(e)
}

#' Print value object and its computational graph
#'
#' @param v A value object
#' @param depth Integer indicates depth (dfs) to recurse down the graph
#' @param indent A character string specifying the indentation for each level in
#' the computational graph (default: "  ")
#'
#' @export
print.value <- function(v, depth = Inf, indent = "  ", ...)
{
  print_helper <- function(x, offset, depth)
  {
    label_str <- if (nchar(x$label) > 0) paste0("'", x$label, "' ") else ""
    dims <- paste0(nrow(x$data), "x", ncol(x$data))
    if (nrow(x$data) == 1 && ncol(x$data) == 1) {
      # Scalar display
      cat(offset, "value(", label_str, "data =", x$data[1,1],
          ", grad =", x$grad[1,1], ")\n", sep = "")
    } else {
      # Matrix display
      cat(offset, "value(", label_str, dims, " matrix)\n", sep = "")
    }
    if (depth > 0)
    {
      for (child in x$prev)
        print_helper(child, paste0(offset, indent), depth - 1)
    }
  }
  print_helper(v, "", depth)
}

#' Check if an object is of class value
#'
#' Determines if the input object is an instance of the value R6 class.
#'
#' @param x The object to be checked
#'
#' @return TRUE if the object is of class value, FALSE otherwise
#' @export
is_value <- function(x) inherits(x, "value")


#' Retrieve the data stored by an object.
#'
#' @param x An object to retrieve the value or data from
#' @param ... additional arguments to pass
#'
#' @return The data matrix of the value object
#' @export
data <- function(x, ...)
{
  UseMethod("data")
}

#' Retrieve the data matrix from a value object
#'
#' @param x A value object
#' @param drop If TRUE (default) and result is 1x1, return scalar.
#'   Set to FALSE to always return a matrix.
#' @param ... additional arguments to pass
#'
#' @return The data as scalar (if 1x1 and drop=TRUE) or matrix
#' @export
data.value <- function(x, drop = TRUE, ...)
{
  if (drop && nrow(x$data) == 1 && ncol(x$data) == 1) {
    return(x$data[1, 1])
  }
  x$data
}

#' Default implementation for retrieving the data from a differentiable object
#'
#' @param x A non-value object
#' @param ... pass additional arguments
#'
#' @return The original object converted to matrix
#' @export
data.default <- function(x, ...)
{
  as_matrix(x)
}


#' Set the data of a value object
#'
#' Assignment function to update the data field of a value object
#' without breaking the computational graph reference.
#'
#' @param x A value object
#' @param value The new numeric value to assign (converted to matrix)
#'
#' @return The modified value object (invisibly)
#'
#' @examples
#' \dontrun{
#' x <- val(5)
#' data(x) <- 10
#' data(x)  # Returns 1x1 matrix containing 10
#' }
#'
#' @export
`data<-` <- function(x, value)
{
  UseMethod("data<-")
}

#' @export
`data<-.value` <- function(x, value)
{
  x$data <- as_matrix(value)
  invisible(x)
}

#' @export
`data<-.default` <- function(x, value)
{
  as_matrix(value)
}


#' Gradient of \code{x} with respect to \code{e} in \code{backward(e)},
#' e.g., dx/de. (applies the chain rule)
#'
#' @param x A differential object
#' @param ... pass additional arguments
#' @return The gradient matrix (same dimensions as data)
#' @export
grad <- function(x, ...)
{
  UseMethod("grad")
}

#' Gradient of a \code{value} object \code{x} with respect to \code{e} in
#' \code{backward(e)}, e.g., dx/de. (applies the chain rule)
#'
#' @param x A value object
#' @param drop If TRUE (default) and result is 1x1, return scalar.
#'   Set to FALSE to always return a matrix.
#' @param ... pass additional arguments
#' @return The gradient as scalar (if 1x1 and drop=TRUE) or matrix
#' @export
grad.value <- function(x, drop = TRUE, ...)
{
  if (drop && nrow(x$grad) == 1 && ncol(x$grad) == 1) {
    return(x$grad[1, 1])
  }
  x$grad
}

#' Default gradient is zero matrix
#'
#' @param x A non-value object
#' @param ... pass additional arguments
#'
#' @return Zero matrix of appropriate dimensions
#' @export
grad.default <- function(x, ...)
{
  m <- as_matrix(x)
  matrix(0, nrow = nrow(m), ncol = ncol(m))
}


#' Default implementation does not propagate gradients. For instance, if we
#' have a constant, then the partial of the constant is not meaningful.
#'
#' @param e A value object for which the backward pass should be performed
#' @param ... pass additional arguments
#'
#' @export
backward.default <- function(x, ...)
{
}

#' @export
zero_grad.default <- function(e)
{

}
