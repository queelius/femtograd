femtograd
================

- <a href="#exponential-distribution"
  id="toc-exponential-distribution">Exponential distribution</a>
- <a href="#automatic-differentiation-ad"
  id="toc-automatic-differentiation-ad">Automatic differentiation (AD)</a>

Here’s a quick demonstration of how to use the
[`femtograd`](https://github.com/queelius/femtograd) R package.

Load it like this:

``` r
library(femtograd)
#> 
#> Attaching package: 'femtograd'
#> The following object is masked from 'package:utils':
#> 
#>     data
```

## Exponential distribution

Let’s create a simple loglikelihood function for the exponential
distribution paramterized by $\lambda$ (failure rate).

We have

$$
  f_{T_i}(t_i | \lambda) = \lambda \exp(-\lambda t_i).
$$

So, the loglikelihood function is just

$$
  \ell(\lambda) = n \log \lambda - \lambda \sum_{i=1}^n t_i.
$$

Let’s generate $n=30$ observations.

``` r
n <- 30
true_rate <- 7.3
data <- rexp(n,true_rate)
head(data)
#> [1] 0.102432677 0.009849689 0.037429092 0.093926112 0.033046038 0.181780460

(mle.rate <- abs(1/mean(data)))
#> [1] 6.245535
```

We see that the MLE $\hat\theta$ is 6.2455349.

# Automatic differentiation (AD)

Finding the value ($\operatorname{argmax}$) that maximizes the
log-likelihood function is trivial to solve in this case, and it has a
closed-form solution. However, to demonstrate the use of `femtograd`, we
will construct a `loglike_exp` function generator that returns an object
that can be automatically differentiated (AD) using backpropogation,
which is an efficient way of applying the chain-rule to expressions
(like $\exp\{y a x^2\}$ using a *computational graph* that represents
the expression.

These kind of computational graphs have the nice property that for any
differentiable expression that we can model in software, its partial
derivative with respect to some node in the graph can be efficiently and
accurately computed without resorting to numerical finite difference
methods or slow, potentially difficult to compose symbolic methods.

There are many libraries that do this. This library itself is based on
the excellent work by Karpathy who developed the Python library known as
[`micrograd`](https://github.com/karpathy/micrograd), which was
developed for the explicit purpose of teaching the basic concept of AD
and backpropagation for minimizing loss functions for neural networks.

Let’s solve for the MLE iteratively as a demonstration of how to use
[`femtograd`](https://github.com/queelius/femtograd). First, we
construct the log-likelihood generator:

``` r
loglike_exp <- function(rate, data)
{
  return(log(rate)*length(data) - rate * sum(data))
}
```

Initially, we guess that $\hat\lambda$ is $1$, which is a terrible
estimate.

``` r
rate <- val(1)
```

Gradient clipping is a technique to prevent taking too large of a step
when gradients become too large (remember that gradients are a *local*
feature, so we generally should not use it to take too big of a step)
during optimization, which can cause instability or overshooting the
optimal value. By limiting the step size, gradient clipping helps ensure
that the optimization takes smaller, more stable steps.

Here is the R code:

``` r
# Takes a gradient `grad` and an optional `max_norm` parameter, which defaults
# to 1. It calculates the gradient's L2 norm (Euclidean norm) and scales the
# gradient down if its norm exceeds the specified max_norm. This is used during
# the gradient ascent loop to help ensure stable optimization.
grad_clip <- function(g, max_norm = 1) {
  norm <- sqrt(sum(g * g))
  if (norm > max_norm) {
    g <- (max_norm / norm) * g
  }
  g
}
```

We find the MLE using a simple iteration (200 loops).

``` r
loglik <- loglike_exp(rate, data)
lr <- 0.2 # learning rate
for (i in 1:200)
{
  zero_grad(loglik)
  backward(loglik)

  data(rate) <- data(rate) + lr * grad_clip(grad(rate))
  if (i %% 50 == 0)
    cat("iteration", i, ", rate =", data(rate), ", drate/dl =", grad(rate), "\n")
}
#> iteration 50 , rate = 6.238781 , drate/dl = 0.006148105 
#> iteration 100 , rate = 6.245533 , drate/dl = 1.447689e-06 
#> iteration 150 , rate = 6.245535 , drate/dl = 3.418377e-10 
#> iteration 200 , rate = 6.245535 , drate/dl = 8.082424e-14
```

Did the gradient ascent method converge to the MLE?

``` r
(converged <- (abs(mle.rate - data(rate)) < 1e-3))
#> [1] TRUE
```

It’s worth pointing out that we did not update `loglik` in the gradient
ascent loop, since we only needed the gradient (score) in this case. If,
however, we had needed to know the log-likelihood for some reason, such
as when using a line search method to avoid overshooting, we would need
to update with `loglik <- loglike_exp(rate, data)` each time through the
loop.
