# Set the data of a value object

Assignment function to update the data field of a value object without
breaking the computational graph reference.

## Usage

``` r
data(x) <- value
```

## Arguments

- x:

  A value object

- value:

  The new numeric value to assign

## Value

The modified value object (invisibly)

## Examples

``` r
if (FALSE) { # \dontrun{
x <- val(5)
data(x) <- 10
data(x)  # Returns 10
} # }
```
