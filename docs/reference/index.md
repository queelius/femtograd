# Package index

## All functions

- [`abs(`*`<value>`*`)`](https://queelius.github.io/femtograd/reference/abs.value.md)
  : Absolute value for value objects

- [`backward()`](https://queelius.github.io/femtograd/reference/backward.md)
  :

  Generic function for the Backward pass for automatic differentiation
  (finds the gradient of every sub-node in the computational graph with
  respect to `e`). In other words, it is responsible for computing the
  gradient with respect to `e`.

- [`backward(`*`<default>`*`)`](https://queelius.github.io/femtograd/reference/backward.default.md)
  : Default implementation does not propagate gradients. For instance,
  if we have a constant, then the partial of the constant is not
  meaningful.

- [`backward(`*`<value>`*`)`](https://queelius.github.io/femtograd/reference/backward.value.md)
  : Backward pass for value objects

- [`confint_mle()`](https://queelius.github.io/femtograd/reference/confint_mle.md)
  : Compute confidence intervals from MLE results

- [`cos(`*`<value>`*`)`](https://queelius.github.io/femtograd/reference/cos.value.md)
  : Cosine function for value objects

- [`` `data<-`() ``](https://queelius.github.io/femtograd/reference/data-set.md)
  : Set the data of a value object

- [`x`](https://queelius.github.io/femtograd/reference/data.md) :
  Retrieve the data stored by an object.

- [`data(`*`<default>`*`)`](https://queelius.github.io/femtograd/reference/data.default.md)
  : Default implementation for retrieving the data from a differentiable
  object

- [`data(`*`<value>`*`)`](https://queelius.github.io/femtograd/reference/data.value.md)
  : Retrieve the value or data from a value object

- [`digamma(`*`<value>`*`)`](https://queelius.github.io/femtograd/reference/digamma.value.md)
  : Digamma (psi) function for value objects

- [`distributions`](https://queelius.github.io/femtograd/reference/distributions.md)
  : Log-likelihood functions for exponential family distributions

- [`` `-`( ``*`<value>`*`)`](https://queelius.github.io/femtograd/reference/dot-value.md)
  : Subtraction for value objects

- [`dual`](https://queelius.github.io/femtograd/reference/dual.md) :
  dual R6 class for forward-mode automatic differentiation

- [`dual_num()`](https://queelius.github.io/femtograd/reference/dual_num.md)
  : Create a dual number

- [`exp(`*`<value>`*`)`](https://queelius.github.io/femtograd/reference/exp.value.md)
  : Exponential function for value objects

- [`find_mle()`](https://queelius.github.io/femtograd/reference/find_mle.md)
  : Find MLE with standard errors

- [`fisher_information()`](https://queelius.github.io/femtograd/reference/fisher_information.md)
  : Compute observed Fisher information matrix

- [`fisher_scoring()`](https://queelius.github.io/femtograd/reference/fisher_scoring.md)
  : Fisher scoring optimizer

- [`grad()`](https://queelius.github.io/femtograd/reference/grad.md) :

  Gradient of `x` with respect to `e` in `backward(e)`, e.g., dx/de.
  (applies the chain rule)

- [`grad(`*`<default>`*`)`](https://queelius.github.io/femtograd/reference/grad.default.md)
  :

  Default gradient is one that does not propograte gradients and is
  zero.`value` object `x` with respect to `e` in `backward(e)`, e.g.,
  dx/de. (applies the chain rule)

- [`grad(`*`<value>`*`)`](https://queelius.github.io/femtograd/reference/grad.value.md)
  :

  Gradient of a `value` object `x` with respect to `e` in `backward(e)`,
  e.g., dx/de. (applies the chain rule)

- [`gradient()`](https://queelius.github.io/femtograd/reference/gradient.md)
  : Compute gradient as a numeric vector

- [`gradient_ascent()`](https://queelius.github.io/femtograd/reference/gradient_ascent.md)
  : Gradient ascent/descent optimizer

- [`gradient_descent()`](https://queelius.github.io/femtograd/reference/gradient_descent.md)
  : Gradient descent (minimize)

- [`hessian()`](https://queelius.github.io/femtograd/reference/hessian.md)
  : Compute Hessian matrix via forward-over-reverse automatic
  differentiation

- [`is_dual()`](https://queelius.github.io/femtograd/reference/is_dual.md)
  : Check if object is a dual number

- [`is_value()`](https://queelius.github.io/femtograd/reference/is_value.md)
  : Check if an object is of class value

- [`lgamma(`*`<value>`*`)`](https://queelius.github.io/femtograd/reference/lgamma.value.md)
  : Log-gamma function for value objects

- [`log(`*`<value>`*`)`](https://queelius.github.io/femtograd/reference/log.value.md)
  : Natural logarithm for value objects

- [`log1p(`*`<value>`*`)`](https://queelius.github.io/femtograd/reference/log1p.value.md)
  : Log(1+x) for value objects

- [`logit()`](https://queelius.github.io/femtograd/reference/logit.md) :
  Logit function for value objects

- [`loglik_bernoulli()`](https://queelius.github.io/femtograd/reference/loglik_bernoulli.md)
  : Bernoulli distribution log-likelihood

- [`loglik_beta()`](https://queelius.github.io/femtograd/reference/loglik_beta.md)
  : Beta distribution log-likelihood

- [`loglik_binomial()`](https://queelius.github.io/femtograd/reference/loglik_binomial.md)
  : Binomial distribution log-likelihood

- [`loglik_exponential()`](https://queelius.github.io/femtograd/reference/loglik_exponential.md)
  : Exponential distribution log-likelihood

- [`loglik_gamma()`](https://queelius.github.io/femtograd/reference/loglik_gamma.md)
  : Gamma distribution log-likelihood

- [`loglik_logistic()`](https://queelius.github.io/femtograd/reference/loglik_logistic.md)
  : Logistic regression log-likelihood (binary)

- [`loglik_negbinom()`](https://queelius.github.io/femtograd/reference/loglik_negbinom.md)
  : Negative binomial log-likelihood

- [`loglik_normal()`](https://queelius.github.io/femtograd/reference/loglik_normal.md)
  : Normal (Gaussian) log-likelihood

- [`loglik_poisson()`](https://queelius.github.io/femtograd/reference/loglik_poisson.md)
  : Poisson distribution log-likelihood

- [`mean(`*`<value>`*`)`](https://queelius.github.io/femtograd/reference/mean.value.md)
  : Mean for value objects

- [`newton_raphson()`](https://queelius.github.io/femtograd/reference/newton_raphson.md)
  : Newton-Raphson optimizer

- [`optimization`](https://queelius.github.io/femtograd/reference/optimization.md)
  : Optimization routines for maximum likelihood estimation

- [`` `+`( ``*`<value>`*`)`](https://queelius.github.io/femtograd/reference/plus-.value.md)
  : Addition for value objects

- [`` `^`( ``*`<value>`*`)`](https://queelius.github.io/femtograd/reference/pow-.value.md)
  : Power operation for value objects.

- [`primal()`](https://queelius.github.io/femtograd/reference/primal.md)
  : Extract primal from dual or return value unchanged

- [`print(`*`<value>`*`)`](https://queelius.github.io/femtograd/reference/print.value.md)
  : Print value object and its computational graph

- [`relu()`](https://queelius.github.io/femtograd/reference/relu.md) :
  ReLU (Rectified Linear Unit) activation function for value objects

- [`sigmoid()`](https://queelius.github.io/femtograd/reference/sigmoid.md)
  : Sigmoid activation function for value objects

- [`sin(`*`<value>`*`)`](https://queelius.github.io/femtograd/reference/sin.value.md)
  : Sine function for value objects

- [`` `/`( ``*`<value>`*`)`](https://queelius.github.io/femtograd/reference/slash-.value.md)
  : Division for value objects

- [`softplus()`](https://queelius.github.io/femtograd/reference/softplus.md)
  : Softplus function for value objects

- [`sqrt(`*`<value>`*`)`](https://queelius.github.io/femtograd/reference/sqrt.value.md)
  : Square root for value objects

- [`std_errors()`](https://queelius.github.io/femtograd/reference/std_errors.md)
  : Compute standard errors from Hessian

- [`sum(`*`<dual>`*`)`](https://queelius.github.io/femtograd/reference/sum.dual.md)
  : Sum for dual numbers

- [`sum(`*`<value>`*`)`](https://queelius.github.io/femtograd/reference/sum.value.md)
  : Summation for value objects

- [`tangent()`](https://queelius.github.io/femtograd/reference/tangent.md)
  : Extract tangent from dual or return 0

- [`tanh(`*`<value>`*`)`](https://queelius.github.io/femtograd/reference/tanh.value.md)
  : Hyperbolic tangent activation function for value objects

- [`` `*`( ``*`<value>`*`)`](https://queelius.github.io/femtograd/reference/times-.value.md)
  : Multiplication for value objects

- [`trigamma(`*`<value>`*`)`](https://queelius.github.io/femtograd/reference/trigamma.value.md)
  : Trigamma function for value objects

- [`val()`](https://queelius.github.io/femtograd/reference/val.md) :

  `value` object constructor

- [`value`](https://queelius.github.io/femtograd/reference/value.md) :
  value R6 class

- [`vcov_matrix()`](https://queelius.github.io/femtograd/reference/vcov_matrix.md)
  : Compute variance-covariance matrix from Hessian

- [`wald_test()`](https://queelius.github.io/femtograd/reference/wald_test.md)
  : Wald test for hypothesis testing
