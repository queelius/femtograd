femtograd
================

Here’s a quick demonstration of how to use the
[`femtograd`](https://github.com/queelius/femtograd) R package.

Load it like this:

``` r
library(femtograd)
```

## Exponential distribution

Let’s create a simple loglikelihood function for the exponential
distribution paramterized by $\lambda$ (failure rate).

We have $$
  f_{T_i}(t_i | \lambda) = \lambda \exp(-\lambda t_i).
$$ So, the loglikelihood function is just $$
  \ell(\lambda) = n \log \lambda - \lambda \sum_{i=1}^n t_i.
$$

Let’s generate $n=30$ observations.

``` r
n <- 30
true_rate <- 7.3
data <- rexp(n,true_rate)
head(data)
#> [1] 0.055221361 0.140126138 0.011210974 0.006967867 0.335711352 0.255494529

(mle.rate <- abs(1/mean(data)))
#> [1] 6.687716
```

We see that the MLE $\hat\theta$ is 6.6877162. Let’s solve for the MLE
iteratively as a demonstration of how to use
[`femtograd`](https://github.com/queelius/femtograd). First, we
construct the log-likelihood function.

``` r
loglike_exp <- function(rate, data)
{
  return(log(rate)*length(data) - rate * sum(data))
}
```

Initially, we guess that $\hat\lambda$ is $5$, which is a terrible
estimate.

``` r
rate <- val(1)
```

Next, we find the MLE using a simple iteration (100 loops).

``` r
# to prevent taking too large of a step
# gradients are a local feature, so we don't want to make too big of jumps
grad_clip <- function(grad, max_norm = 1)
{
  norm <- sqrt(sum(grad * grad))
  if (norm > max_norm)
    grad <- (max_norm / norm) * grad
  grad
}

loglik <- loglike_exp(rate, data)
lr <- 0.2
for (i in 1:200)
{
  zero_grad(loglik)
  backward(loglik)

  rate$data <- rate$data + lr * grad_clip(rate$grad)
  if (i %% 10 == 0)
    cat("iteration", i, ", rate =", rate$data, ", rate.grad =", rate$grad, "\n")
}
#> iteration 10 , rate = 3 , rate.grad = 6.228449 
#> iteration 20 , rate = 5 , rate.grad = 1.764164 
#> iteration 30 , rate = 6.338889 , rate.grad = 0.2906576 
#> iteration 40 , rate = 6.60888 , rate.grad = 0.06205048 
#> iteration 50 , rate = 6.66924 , rate.grad = 0.01436616 
#> iteration 60 , rate = 6.683351 , rate.grad = 0.003384349 
#> iteration 70 , rate = 6.686683 , rate.grad = 0.0008004893 
#> iteration 80 , rate = 6.687472 , rate.grad = 0.0001895165 
#> iteration 90 , rate = 6.687658 , rate.grad = 4.487825e-05 
#> iteration 100 , rate = 6.687703 , rate.grad = 1.062791e-05 
#> iteration 110 , rate = 6.687713 , rate.grad = 2.516895e-06 
#> iteration 120 , rate = 6.687715 , rate.grad = 5.960515e-07 
#> iteration 130 , rate = 6.687716 , rate.grad = 1.411571e-07 
#> iteration 140 , rate = 6.687716 , rate.grad = 3.342887e-08 
#> iteration 150 , rate = 6.687716 , rate.grad = 7.916637e-09 
#> iteration 160 , rate = 6.687716 , rate.grad = 1.874821e-09 
#> iteration 170 , rate = 6.687716 , rate.grad = 4.439951e-10 
#> iteration 180 , rate = 6.687716 , rate.grad = 1.05147e-10 
#> iteration 190 , rate = 6.687716 , rate.grad = 2.490097e-11 
#> iteration 200 , rate = 6.687716 , rate.grad = 5.896617e-12
```

Did it converge to the MLE?

``` r
(converged <- (abs(mle.rate - rate$data) < 1e-3))
#> [1] TRUE
```

It’s worth pointing out that we did not update `loglik` in the loop,
since we only needed the gradient (score) to do that particular gradient
ascent. If, however, we had used the log-likelihood value at `rate`, for
instance using a line search method to avoid having to specify a step
size, then we would need to update `loglik` each time through the loop.
