#' Log-likelihood functions for exponential family distributions
#'
#' These functions compute log-likelihoods that can be differentiated
#' using femtograd's automatic differentiation. Each function returns
#' a value object suitable for gradient-based optimization.
#'
#' @name distributions
NULL


#' Normal (Gaussian) log-likelihood
#'
#' Computes the log-likelihood for i.i.d. normal observations.
#' L(μ,σ²|x) = -n/2 log(2π) - n/2 log(σ²) - Σ(xᵢ-μ)²/(2σ²)
#'
#' @param mu Mean parameter (value object)
#' @param sigma Standard deviation parameter (value object), must be positive
#' @param x Numeric vector of observations
#'
#' @return A value object representing the log-likelihood
#'
#' @examples
#' \dontrun{
#' x <- rnorm(100, mean = 5, sd = 2)
#' mu <- val(0)
#' sigma <- val(1)
#' ll <- loglik_normal(mu, sigma, x)
#' backward(ll)
#' grad(mu)  # score for mu
#' }
#'
#' @export
loglik_normal <- function(mu, sigma, x)
{
  n <- length(x)
  sum_x <- sum(x)
  sum_x2 <- sum(x^2)

  # -n/2 * log(2*pi) - n*log(sigma) - sum((x-mu)^2) / (2*sigma^2)
  # = -n/2 * log(2*pi) - n*log(sigma) - (sum_x2 - 2*mu*sum_x + n*mu^2) / (2*sigma^2)
  const <- -n / 2 * log(2 * pi)
  sigma2 <- sigma * sigma

  # Build expression starting from parameters to ensure correct dispatch
  # log_term = -n * log(sigma) = log(sigma) * (-n)
  log_term <- log(sigma) * (-n)

  # quadratic_term = (sum_x2 - 2*mu*sum_x + n*mu^2) / (2*sigma^2)
  quadratic <- (mu * mu * n - mu * (2 * sum_x) + sum_x2) / (sigma2 * 2)

  # Return: const + log_term - quadratic
  # Use log_term + const to ensure dispatch on value/dual first
  log_term + const - quadratic
}


#' Exponential distribution log-likelihood
#'
#' Computes the log-likelihood for i.i.d. exponential observations.
#' L(λ|x) = n*log(λ) - λ*Σxᵢ
#'
#' @param rate Rate parameter λ (value object), must be positive
#' @param x Numeric vector of observations (must be non-negative)
#'
#' @return A value object representing the log-likelihood
#'
#' @examples
#' \dontrun{
#' x <- rexp(50, rate = 2)
#' rate <- val(1)
#' ll <- loglik_exponential(rate, x)
#' backward(ll)
#' grad(rate)  # should be n/rate - sum(x)
#' }
#'
#' @export
loglik_exponential <- function(rate, x)
{
  n <- length(x)
  sum_x <- sum(x)
  log(rate) * n - rate * sum_x
}


#' Poisson distribution log-likelihood
#'
#' Computes the log-likelihood for i.i.d. Poisson observations.
#' L(λ|x) = Σxᵢ*log(λ) - n*λ - Σlog(xᵢ!)
#'
#' @param lambda Rate parameter λ (value object), must be positive
#' @param x Integer vector of observations (counts)
#'
#' @return A value object representing the log-likelihood
#'
#' @details
#' The term Σlog(xᵢ!) is constant w.r.t. λ and is included for completeness.
#'
#' @examples
#' \dontrun{
#' x <- rpois(100, lambda = 3)
#' lambda <- val(1)
#' ll <- loglik_poisson(lambda, x)
#' backward(ll)
#' # MLE is mean(x)
#' }
#'
#' @export
loglik_poisson <- function(lambda, x)
{
  n <- length(x)
  sum_x <- sum(x)
  sum_logfact <- sum(lgamma(x + 1))  # log(x!) = lgamma(x+1)

  log(lambda) * sum_x - lambda * n - sum_logfact
}


#' Binomial distribution log-likelihood
#'
#' Computes the log-likelihood for binomial observations.
#' L(p|x,n) = Σ[xᵢ*log(p) + (nᵢ-xᵢ)*log(1-p) + log(C(nᵢ,xᵢ))]
#'
#' @param p Success probability (value object), must be in (0,1)
#' @param x Integer vector of successes
#' @param size Integer vector of trial counts (or single value if constant)
#'
#' @return A value object representing the log-likelihood
#'
#' @examples
#' \dontrun{
#' x <- rbinom(50, size = 10, prob = 0.3)
#' p <- val(0.5)
#' ll <- loglik_binomial(p, x, size = 10)
#' backward(ll)
#' # MLE is sum(x) / sum(size)
#' }
#'
#' @export
loglik_binomial <- function(p, x, size)
{
  if (length(size) == 1) size <- rep(size, length(x))

  sum_x <- sum(x)
  sum_n_minus_x <- sum(size - x)
  const <- sum(lchoose(size, x))  # binomial coefficients (constant)

  log(p) * sum_x + log(val(1) - p) * sum_n_minus_x + const
}


#' Bernoulli distribution log-likelihood
#'
#' Special case of binomial with size=1.
#' L(p|x) = Σ[xᵢ*log(p) + (1-xᵢ)*log(1-p)]
#'
#' @param p Success probability (value object), must be in (0,1)
#' @param x Binary vector (0 or 1)
#'
#' @return A value object representing the log-likelihood
#'
#' @export
loglik_bernoulli <- function(p, x)
{
  sum_x <- sum(x)
  n <- length(x)
  log(p) * sum_x + log(val(1) - p) * (n - sum_x)
}


#' Gamma distribution log-likelihood
#'
#' Computes the log-likelihood for i.i.d. gamma observations.
#' L(α,β|x) = n*α*log(β) - n*log(Γ(α)) + (α-1)*Σlog(xᵢ) - β*Σxᵢ
#'
#' @param shape Shape parameter α (value object), must be positive
#' @param rate Rate parameter β (value object), must be positive
#' @param x Numeric vector of observations (must be positive)
#'
#' @return A value object representing the log-likelihood
#'
#' @details
#' The gamma distribution is parameterized with shape α and rate β,
#' where E[X] = α/β and Var[X] = α/β².
#'
#' @export
loglik_gamma <- function(shape, rate, x)
{
  n <- length(x)
  sum_x <- sum(x)
  sum_log_x <- sum(log(x))

  shape * log(rate) * n - lgamma(shape) * n +
    (shape - 1) * sum_log_x - rate * sum_x
}


#' Beta distribution log-likelihood
#'
#' Computes the log-likelihood for i.i.d. beta observations.
#' L(α,β|x) = n*[log(Γ(α+β)) - log(Γ(α)) - log(Γ(β))] +
#'            (α-1)*Σlog(xᵢ) + (β-1)*Σlog(1-xᵢ)
#'
#' @param alpha Shape parameter α (value object), must be positive
#' @param beta Shape parameter β (value object), must be positive
#' @param x Numeric vector of observations in (0,1)
#'
#' @return A value object representing the log-likelihood
#'
#' @export
loglik_beta <- function(alpha, beta, x)
{
  n <- length(x)
  sum_log_x <- sum(log(x))
  sum_log_1mx <- sum(log(1 - x))

  (lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta)) * n +
    (alpha - 1) * sum_log_x + (beta - 1) * sum_log_1mx
}


#' Negative binomial log-likelihood
#'
#' Computes the log-likelihood for negative binomial (failures until r successes).
#' L(r,p|x) = Σ[log(C(xᵢ+r-1, xᵢ)) + r*log(p) + xᵢ*log(1-p)]
#'
#' @param r Number of successes parameter (value object or fixed positive)
#' @param p Success probability (value object), must be in (0,1)
#' @param x Integer vector of failure counts
#'
#' @return A value object representing the log-likelihood
#'
#' @export
loglik_negbinom <- function(r, p, x)
{
  n <- length(x)
  sum_x <- sum(x)

  # Binomial coefficient term: log(C(x+r-1, x)) = lgamma(x+r) - lgamma(x+1) - lgamma(r)
  if (is_value(r)) {
    binom_term <- sum(lgamma(x + data(r))) - sum(lgamma(x + 1)) - lgamma(r) * n
  } else {
    binom_term <- sum(lgamma(x + r)) - sum(lgamma(x + 1)) - lgamma(r) * n
  }

  binom_term + r * log(p) * n + log(1 - p) * sum_x
}


#' Logistic regression log-likelihood (binary)
#'
#' Computes the log-likelihood for binary logistic regression.
#' L(β|X,y) = Σ[yᵢ*log(pᵢ) + (1-yᵢ)*log(1-pᵢ)]
#' where pᵢ = sigmoid(Xᵢ·β)
#'
#' @param beta Coefficient vector (list of value objects)
#' @param X Design matrix (n x p numeric matrix)
#' @param y Binary response vector (0 or 1)
#'
#' @return A value object representing the log-likelihood
#'
#' @details
#' Uses the numerically stable form:
#' log(p) = -log(1 + exp(-η)) and log(1-p) = -log(1 + exp(η))
#' where η = Xβ
#'
#' @export
loglik_logistic <- function(beta, X, y)
{
  # Compute linear predictor η = Xβ
  n <- nrow(X)
  p <- length(beta)

  # Initialize eta as val(0) for each observation
  eta_list <- lapply(1:n, function(i) {
    # η_i = Σ_j X[i,j] * β_j
    terms <- lapply(1:p, function(j) X[i, j] * beta[[j]])
    Reduce(`+`, terms)
  })

  # Log-likelihood: Σ[y*η - log(1 + exp(η))]
  # This is numerically stable form
  ll_terms <- mapply(function(eta, yi) {
    yi * eta - log(1 + exp(eta))
  }, eta_list, y, SIMPLIFY = FALSE)

  Reduce(`+`, ll_terms)
}


#' Weibull distribution log-likelihood
#'
#' Computes the log-likelihood for i.i.d. Weibull observations.
#' L(k, λ | x) = n*log(k) - n*k*log(λ) + (k-1)*Σlog(xᵢ) - Σ(xᵢ/λ)^k
#'
#' @param shape Shape parameter k (value object), must be positive
#' @param scale Scale parameter λ (value object), must be positive
#' @param x Numeric vector of observations (must be positive)
#'
#' @return A value object representing the log-likelihood
#'
#' @details
#' The Weibull distribution is commonly used in survival analysis and
#' reliability engineering. It generalizes the exponential distribution
#' (k=1 gives exponential with rate 1/λ).
#'
#' @examples
#' \dontrun{
#' x <- rweibull(100, shape = 2, scale = 3)
#'
#' # Using log-parameterization for positivity
#' result <- fit(
#'   function(log_shape, log_scale) {
#'     shape <- exp(log_shape)
#'     scale <- exp(log_scale)
#'     loglik_weibull(shape, scale, x)
#'   },
#'   params = c(log_shape = 0, log_scale = 0)
#' )
#' }
#'
#' @export
loglik_weibull <- function(shape, scale, x) {
  n <- length(x)
  sum_log_x <- sum(log(x))

  # L = n*log(k) - n*k*log(λ) + (k-1)*Σlog(xᵢ) - Σ(xᵢ/λ)^k
  # For autodiff, we need to build the expression from parameters

  # First term: n * log(shape)
  term1 <- log(shape) * n

  # Second term: -n * shape * log(scale)
  term2 <- shape * log(scale) * (-n)

  # Third term: (shape - 1) * sum(log(x))
  # Need to handle (shape - 1) carefully for autodiff
  term3 <- (shape - 1) * sum_log_x

  # Fourth term: -sum((x/scale)^shape)
  # This is tricky for autodiff. We need:
  # -sum(x^shape) / scale^shape = -sum(x^shape) * scale^(-shape)
  # For stability, use log-space: exp(shape * log(x) - shape * log(scale))

  # Compute -Σ(xᵢ/λ)^k = -Σexp(k*log(xᵢ) - k*log(λ))
  # = -exp(log(Σexp(k*log(xᵢ) - k*log(λ))))
  # But sum of exponentials doesn't simplify nicely.

  # Alternative: compute as scalar multiplication
  # -Σ(xᵢ/λ)^k where shape is a value object
  # If shape is constant, we could precompute x^shape
  # But for full autodiff, we need to handle this differently

  # For practical cases, if scale is the main parameter:
  if (is_value(shape) || is_dual(shape)) {
    # Use differentiable form
    # Compute sum of (x/scale)^shape = sum(exp(shape * log(x/scale)))
    # = sum(exp(shape * (log(x) - log(scale))))
    # = sum(exp(shape * log(x) - shape * log(scale)))
    # We sum over i: exp(shape * log(x[i]) - shape * log(scale))

    # This is expensive but necessary for full differentiability
    log_x <- log(x)
    # term4 = -sum_i exp(shape * log(x[i]) - shape * log(scale))
    #       = -sum_i exp(shape * (log(x[i]) - log(scale)))

    # Build as a sum over observations
    exp_terms <- lapply(log_x, function(lxi) {
      exp(shape * (lxi - log(scale)))
    })
    term4_sum <- Reduce(`+`, exp_terms)
    term4 <- -term4_sum
  } else {
    # shape is numeric, can precompute
    x_pow_k <- sum(x^shape)
    term4 <- -(x_pow_k) * scale^(-shape)
  }

  term1 + term2 + term3 + term4
}


#' Pareto distribution log-likelihood
#'
#' Computes the log-likelihood for i.i.d. Pareto observations.
#' L(α, xₘ | x) = n*log(α) + n*α*log(xₘ) - (α+1)*Σlog(xᵢ)
#'
#' @param alpha Shape parameter α (value object), must be positive
#' @param x_min Minimum/scale parameter xₘ (fixed positive number).
#'   All observations must be >= x_min.
#' @param x Numeric vector of observations (must be >= x_min)
#'
#' @return A value object representing the log-likelihood
#'
#' @details
#' The Pareto distribution is used to model heavy-tailed phenomena like
#' income distributions, city sizes, etc. Here x_min is typically known
#' (e.g., min(x)) and alpha is estimated.
#'
#' @examples
#' \dontrun{
#' # Generate Pareto data
#' alpha_true <- 2
#' x_min <- 1
#' u <- runif(100)
#' x <- x_min * (1 - u)^(-1/alpha_true)
#'
#' # Fit (alpha only, x_min = min(x) is fixed)
#' result <- fit(
#'   function(log_alpha) {
#'     alpha <- exp(log_alpha)
#'     loglik_pareto(alpha, x_min = min(x), x)
#'   },
#'   params = c(log_alpha = 0)
#' )
#' }
#'
#' @export
loglik_pareto <- function(alpha, x_min, x) {
  n <- length(x)
  sum_log_x <- sum(log(x))

  # L = n*log(α) + n*α*log(xₘ) - (α+1)*Σlog(xᵢ)
  log(alpha) * n + alpha * (log(x_min) * n) - (alpha + 1) * sum_log_x
}
