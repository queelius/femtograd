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
distribution paramterized by
![\\lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda
"\\lambda") (failure rate).

We have   
![&#10; f\_{T\_i}(t\_i | \\lambda) = \\lambda \\exp(-\\lambda
t\_i).&#10;](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%20%20f_%7BT_i%7D%28t_i%20%7C%20%5Clambda%29%20%3D%20%5Clambda%20%5Cexp%28-%5Clambda%20t_i%29.%0A
"
  f_{T_i}(t_i | \\lambda) = \\lambda \\exp(-\\lambda t_i).
")  
So, the loglikelihood function is just   
![&#10; \\ell(\\lambda) = n \\log \\lambda - \\lambda \\sum\_{i=1}^n
t\_i.&#10;](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%20%20%5Cell%28%5Clambda%29%20%3D%20n%20%5Clog%20%5Clambda%20-%20%5Clambda%20%5Csum_%7Bi%3D1%7D%5En%20t_i.%0A
"
  \\ell(\\lambda) = n \\log \\lambda - \\lambda \\sum_{i=1}^n t_i.
")  

Let’s generate
![n=30](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%3D30
"n=30") observations.

``` r
n <- 30
true_rate <- 7.3
data <- rexp(n,true_rate)
head(data)
#> [1] 0.05387135 0.07232008 0.12197395 0.03077618 0.17593061 0.09581282

(mle.rate <- abs(1/mean(data)))
#> [1] 8.687603
```

We see that the MLE
![\\hat\\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%5Ctheta
"\\hat\\theta") is 8.6876029. Let’s solve for the MLE iteratively as a
demonstration of how to use
[`femtograd`](https://github.com/queelius/femtograd). First, we
construct the log-likelihood function.

``` r
loglike_exp <- function(rate, data)
{
  return(log(rate)*length(data) - rate * sum(data))
}
```

Initially, we guess that
![\\hat\\lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%5Clambda
"\\hat\\lambda") is
![5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;5
"5"), which is a terrible estimate.

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
#> iteration 10 , rate = 3 , rate.grad = 7.261089 
#> iteration 20 , rate = 5 , rate.grad = 2.796803 
#> iteration 30 , rate = 6.991714 , rate.grad = 0.9585682 
#> iteration 40 , rate = 8.038382 , rate.grad = 0.3076895 
#> iteration 50 , rate = 8.417195 , rate.grad = 0.1212328 
#> iteration 60 , rate = 8.57176 , rate.grad = 0.05082324 
#> iteration 70 , rate = 8.63742 , rate.grad = 0.02181832 
#> iteration 80 , rate = 8.665762 , rate.grad = 0.009459143 
#> iteration 90 , rate = 8.678078 , rate.grad = 0.004118187 
#> iteration 100 , rate = 8.683446 , rate.grad = 0.001796177 
#> iteration 110 , rate = 8.685788 , rate.grad = 0.0007840346 
#> iteration 120 , rate = 8.68681 , rate.grad = 0.0003423505 
#> iteration 130 , rate = 8.687257 , rate.grad = 0.0001495106 
#> iteration 140 , rate = 8.687452 , rate.grad = 6.529823e-05 
#> iteration 150 , rate = 8.687537 , rate.grad = 2.851959e-05 
#> iteration 160 , rate = 8.687574 , rate.grad = 1.245635e-05 
#> iteration 170 , rate = 8.68759 , rate.grad = 5.44052e-06 
#> iteration 180 , rate = 8.687597 , rate.grad = 2.376245e-06 
#> iteration 190 , rate = 8.6876 , rate.grad = 1.037869e-06 
#> iteration 200 , rate = 8.687602 , rate.grad = 4.533085e-07
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
