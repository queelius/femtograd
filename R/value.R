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
#' @field op A character string representing the operation associated with the node
#'
#' @import R6
#' @export
value <- R6Class("value",
   public = list(
     data = NULL,
     grad = 0,
     backward_fn = NULL,
     prev = NULL,
     op = "",

     #' @description Initializes a new value object with the given data, list of children, and operation name.
     #' @param data Numeric value of the object
     #' @param children List of previous nodes in the computational graph (default: empty list)
     #' @param op A character string representing the operation associated with the node (default: empty string)
     initialize = function(data, children = list(), op = "")
     {
       self$data <- data
       self$grad <- 0
       self$backward_fn <- function() NULL
       self$prev <- children
       self$op <- op
     })
)

#' value constructor
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

#' Representation of a value object
#'
#' Creates a string representation of a value object including its data and gradient values.
#'
#' @param x A value object
#'
#' @return A string representation of the value object
#' @export
repr <- function(x) sprintf("value(data=%s, grad=%s)", x$data, x$grad)

#' Backward pass for automatic differentiation
#'
#' Generic function for the backward pass, responsible for computing gradients in the
#' computational graph. This function should be implemented for specific classes.
#'
#' @param e An object for which the backward pass should be performed
#'
#' @export
backward <- function(e) UseMethod("backward")

#' Backward pass for value objects
#'
#' Performs the backward pass (gradient computation) for a value object in the
#' computational graph with automatic differentiation.
#'
#' @param e A value object for which the backward pass should be performed
#'
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
#' @param indent A character string specifying the indentation for each level in
#' the computational graph (default: "  ")
#'
#' @export
print.value <- function(v, indent = "  ")
{
  print_tree <- function(x, offset)
  {
    cat(offset, "value( data =", x$data, ", grad =", x$grad, ", op =", x$op, ")\n")
    for (child in x$prev)
      print_tree(child, paste0(offset, indent))
  }
  print_tree(v, "")
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


#' Retrieve the value or data from an object
#'
#' @param x An object to retrieve the value or data from
#'
#' @return The value or data of the object if it is a value object, or the original object if it is not a value object
#' @export
retrieve <- function(x)
{
  UseMethod("retrieve")
}

#' Retrieve the value or data from a value object
#'
#' @param x A value object
#'
#' @return The value or data of the value object
#' @export
retrieve.value <- function(x)
{
  x$data
}

#' Retrieve the value or data from a non-value object
#'
#' @param x A non-value object
#'
#' @return The original object
#' @export
retrieve.default <- function(x)
{
  x
}
