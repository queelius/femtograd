femtograd
================

  - [0.1 Exponential distribution](#01-exponential-distribution)

Here’s a quick demonstration of how to use the
[`femtograd`](https://github.com/queelius/femtograd) R package.

Load it like this:

``` r
library(femtograd)
```

## 0.1 Exponential distribution

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
#> [1] 0.067121114 0.022125584 0.141639518 0.001758898 0.025229679 0.138780163

(mle.rate <- abs(1/mean(data)))
#> [1] 10.50944
```

We see that the MLE
![\\hat\\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%5Ctheta
"\\hat\\theta") is 10.5094428. Let’s solve for the MLE iteratively as a
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
#> iteration 10 , rate = 3 , rate.grad = 7.85971 
#> iteration 20 , rate = 5 , rate.grad = 3.395424 
#> iteration 30 , rate = 7 , rate.grad = 1.557189 
#> iteration 40 , rate = 8.749491 , rate.grad = 0.6237931 
#> iteration 50 , rate = 9.584124 , rate.grad = 0.2949889 
#> iteration 60 , rate = 10.00197 , rate.grad = 0.154106 
#> iteration 70 , rate = 10.22561 , rate.grad = 0.08406503 
#> iteration 80 , rate = 10.34909 , rate.grad = 0.04685764 
#> iteration 90 , rate = 10.41835 , rate.grad = 0.02641968 
#> iteration 100 , rate = 10.45754 , rate.grad = 0.01499036 
#> iteration 110 , rate = 10.47982 , rate.grad = 0.008535471 
#> iteration 120 , rate = 10.49252 , rate.grad = 0.004869762 
#> iteration 130 , rate = 10.49977 , rate.grad = 0.002781499 
#> iteration 140 , rate = 10.50391 , rate.grad = 0.001589754 
#> iteration 150 , rate = 10.50628 , rate.grad = 0.0009089512 
#> iteration 160 , rate = 10.50763 , rate.grad = 0.0005198073 
#> iteration 170 , rate = 10.50841 , rate.grad = 0.000297301 
#> iteration 180 , rate = 10.50885 , rate.grad = 0.0001700514 
#> iteration 190 , rate = 10.5091 , rate.grad = 9.727046e-05 
#> iteration 200 , rate = 10.50925 , rate.grad = 5.564057e-05
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
