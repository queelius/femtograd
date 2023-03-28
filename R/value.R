#' value R6 class
#'
#' Represents a node in the computational graph with automatic differentiation.
#' Each value object holds data, gradient, a backward function, previous nodes,
#' and an operation name.
#'
#' @field data Numeric value of the object
#' @field grad Gradient value, initially set to 0
#' @field backward_fn A function that performs the backward pass (gradient computation)
#' @field prev A list of previous nodes in the computational graph
#'
#' @import R6
#' @export
value <- R6Class("value",
   public = list(
     data = NULL,
     grad = 0,
     backward_fn = NULL,
     prev = NULL,

     #' @description Initializes a new value object with the given data, list of children, and operation name.
     #' @param data Numeric value of the object
     #' @param children List of previous nodes in the computational graph (default: empty list)
     initialize = function(data, children = list())
     {
       self$data <- data
       self$grad <- 0
       self$backward_fn <- function() NULL
       self$prev <- children
     })
)

#' \code{value} object constructor
#'
#' @param data Numeric value of the object
#'
#' Represents a value node in the computational graph with automatic
#' differentiation. This constructor function creates a new value object with
#' the specified data.
#' @export
val <- function(data)
{
  value$new(data)
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
#' computational graph with automatic differentiation.
#'
#' @param e A value object for which the backward pass should be performed
#' @param ... pass additional arguments
#' @export
backward.value <- function(x)
{
  topo <- topological_sort(x)
  x$grad <- 1
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
    v$grad <- 0
    for (child in v$prev)
      reset_grad(child)
  }
  reset_grad(e)
}

#' Print value object and its computational graph
#'
#' @param v A value object
#' @param depth Integer indicates depth (dfs) to recursve down the graph
#' @param indent A character string specifying the indentation for each level in
#' the computational graph (default: "  ")
#'
#' @export
print.value <- function(v, depth = Inf, indent = "  ")
{
  print_helper <- function(x, offset, depth)
  {
    cat(offset, "value( data =", x$data, ", grad =", x$grad, ")\n")
    if (depth > 0)
    {
      for (child in x$prev)
        print_helper(child, depth-1, paste0(offset, indent))
    }
  }
  print_helper(v, depth, "")
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
#' @return The value or data of a differnetiable object
#' @export
data <- function(x, ...)
{
  UseMethod("data")
}

#' Retrieve the value or data from a value object
#'
#' @param x A value object
#' @param ... additional arguments to pass
#'
#' @return The value or data of the value object
#' @export
data.value <- function(x, ...)
{
  x$data
}

#' Default implementation for retrieving the data from a differentiable object
#'
#' @param x A non-value object
#' @param ... pass additional arguments
#'
#' @return The original object
#' @export
data.default <- function(x, ...)
{
  x
}


#' Gradient of \code{x} with respect to \code{e} in \code{backward(e)},
#' e.g., dx/de. (applies the chain rule)
#'
#' @param x A differential object
#' @param ... pass additional arguments
#' @return The gradient as previously defined
#' @export
grad <- function(x, ...)
{
  UseMethod("grad")
}

#' Gradient of a \code{value} object \code{x} with respect to \code{e} in
#' \code{backward(e)}, e.g., dx/de. (applies the chain rule)
#'
#' @param x A value object
#' @param ... pass additional arguments
#' @return The value or data of the value object
#' @export
grad.value <- function(x, ...)
{
  x$grad
}

#' Default gradient is one that does not propograte gradients and is zero.\code{value} object \code{x} with respect to \code{e} in
#' \code{backward(e)}, e.g., dx/de. (applies the chain rule)
#'
#' @param x A value object
#' @param ... pass additional arguments
#'
#' @return The value or data of the value object
#' @export
grad.default <- function(x, ...)
{
  0.0
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
