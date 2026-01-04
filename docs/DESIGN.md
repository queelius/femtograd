# femtograd Design Decisions

## Matrix-Based Representation

All values in femtograd use matrix representation internally. This
document explains why.

### The Problem

Early versions of femtograd had inconsistent return types:

``` r
x <- val(3)
y <- val(c(1, 2, 3))
z <- x + c(1, 2, 3)  # scalar val + vector

data(x)    # scalar
data(y)    # vector
data(z)    # vector (from broadcasting)
grad(x)    # scalar or vector, depending on computation graph
```

This unpredictability made the API difficult to use correctly. Users
couldn’t know whether
[`grad()`](https://queelius.github.io/femtograd/reference/grad.md) would
return a scalar, vector, or something else without tracing through the
entire computation.

### The Solution

Standardize on matrices. Everything is 2-dimensional:

- Scalars are 1×1 matrices
- Column vectors are n×1 matrices
- Row vectors are 1×n matrices
- Matrices are m×n matrices

This means: -
[`data()`](https://queelius.github.io/femtograd/reference/data.md)
always returns a matrix -
[`grad()`](https://queelius.github.io/femtograd/reference/grad.md)
always returns a matrix (same dimensions as the value) -
[`hessian()`](https://queelius.github.io/femtograd/reference/hessian.md)
always returns a matrix (p×p for p parameters)

### Why Matrices Instead of Tensors?

We considered full tensor support (arbitrary n-dimensional arrays) but
rejected it:

**1. Statistical applications don’t need it**

The target use cases are maximum likelihood estimation, Fisher
information, Hessian computation, and gradient-based optimization. These
involve:

- Scalar loss functions
- Vector parameters (p×1)
- Hessian matrices (p×p)
- Variance-covariance matrices (p×p)

No standard statistical computation requires 3+ dimensional tensors.

**2. Pedagogy over generality**

femtograd exists to teach automatic differentiation concepts clearly.
Tensors introduce:

- Shape algebra and broadcasting rules
- Stride calculations and memory layout
- View vs. copy semantics
- Einsum notation for contractions

These obscure the core AD concepts (computational graphs, chain rule,
forward/reverse modes) behind implementation complexity.

**3. R already has tensor autodiff**

The `torch` package provides GPU-accelerated tensors with full autodiff
support. Users needing tensor operations should use it. femtograd is not
competing with PyTorch—it’s a teaching tool for statistics.

**4. Implementation simplicity**

Matrix operations in R are well-supported and fast (BLAS/LAPACK).
Implementing a custom tensor library would be substantial work for
little benefit in our use case.

### Dimension Conventions

We follow standard mathematical conventions:

- **Parameters**: Column vectors (p×1)
- **Gradients**: Same shape as the value (p×1 for parameter vectors)
- **Jacobians**: If f: ℝⁿ → ℝᵐ, then J is m×n
- **Hessians**: Square matrices (p×p)
- **Fisher information**: Square matrices (p×p)

### Implications for the API

``` r
# Single parameter (1×1)
theta <- val(matrix(5, 1, 1))

# Multiple parameters (p×1 column vector)
params <- val(matrix(c(mu, sigma), ncol = 1))

# Gradient has same shape as value
grad(params)  # p×1 matrix

# Hessian is p×p
hessian(loss_fn, params)  # p×p matrix
```

### Trade-offs

**What we gain:** - Predictable API behavior - Clear dimension
semantics - Alignment with mathematical notation - Simpler
implementation

**What we lose:** - Convenience of implicit scalar handling - Some
R-idiomatic flexibility (e.g., `val(5)` without explicit matrix) -
Higher-dimensional operations (rarely needed in statistics)

### Compatibility Notes

For convenience,
[`val()`](https://queelius.github.io/femtograd/reference/val.md) accepts
scalars and vectors and converts them to matrices internally:

``` r
val(5)           # becomes 1×1 matrix
val(c(1, 2, 3))  # becomes 3×1 column vector
```

This preserves a simple API while maintaining internal consistency.

**Smart dropping**: By default,
[`data()`](https://queelius.github.io/femtograd/reference/data.md) and
[`grad()`](https://queelius.github.io/femtograd/reference/grad.md)
return scalars for 1×1 matrices (using `drop=TRUE`). This means most
scalar operations work naturally:

``` r
x <- val(3)
y <- x^2
backward(y)
data(y)   # Returns 9 (scalar), not matrix(9, 1, 1)
grad(x)   # Returns 6 (scalar)
```

For explicit matrix returns, use `drop=FALSE`:

``` r
data(x, drop = FALSE)  # Returns 1×1 matrix
```

### Broadcasting

Scalar (1×1) matrices broadcast over larger matrices in element-wise
operations:

``` r
a <- val(2)           # 1×1 matrix
x <- val(c(1, 2, 3))  # 3×1 matrix
y <- a * x            # 3×1 matrix: c(2, 4, 6)

backward(y)
grad(a)  # Scalar: sum of all gradients = 6
grad(x)  # 3×1 matrix: c(2, 2, 2)
```

This is essential for statistical operations like `lambda * x` in
exponential distributions.
