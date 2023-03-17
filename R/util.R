#' Set membership test for value objects.
#'
#' @param xs A list of value objects
#' @param obj A value object
#' @return Boolean, determine if \code{xs} has \code{obj}
#' @keywords internal
contains <- function(xs, obj)
{
  any(sapply(xs, identical, obj))
}

#' Build the topological order for a given value object
#'
#' @param v A value object
#' @return A list representing the topological order of nodes
#' @keywords internal
topological_sort <- function(v)
{
  res <- list()
  visited <- list()

  build_topo <- function(x)
  {
    if (!contains(visited, x))
    {
      visited <<- append(visited, x)
      for (child in x$prev)
        build_topo(child)
      res <<- append(res, x)
    }
  }
  build_topo(v)
  res
}
