#' Hypothesis Testing for Fitted Models
#'
#' Provides hypothesis testing functions that produce results compatible
#' with the hypothesize package interface. Results have `stat`, `p.value`,
#' and `dof` components accessible via standard methods.
#'
#' @name hypothesis_tests
NULL


#' Likelihood Ratio Test
#'
#' Computes the likelihood ratio test (LRT) for comparing nested models.
#' Works with femtofit objects or raw log-likelihood values.
#'
#' @param null_model Either a `femtofit` object (simpler model) or a numeric
#'   log-likelihood value.
#' @param alt_model Either a `femtofit` object (more complex model) or a numeric
#'   log-likelihood value.
#' @param df Degrees of freedom. If NULL and both arguments are femtofit objects,
#'   computed as the difference in number of parameters.
#'
#' @return A `likelihood_ratio_test` object (also inherits from `hypothesis_test`)
#'   containing:
#'   \describe{
#'     \item{stat}{The LRT statistic: -2 * (loglik_null - loglik_alt)}
#'     \item{p.value}{P-value from chi-squared distribution}
#'     \item{dof}{Degrees of freedom}
#'     \item{null_loglik}{Log-likelihood of null model}
#'     \item{alt_loglik}{Log-likelihood of alternative model}
#'   }
#'
#' @details
#' The likelihood ratio test statistic is:
#' \deqn{\Lambda = -2 (\ell_0 - \ell_1)}
#' where \eqn{\ell_0} is the log-likelihood under the null (simpler) model
#' and \eqn{\ell_1} is the log-likelihood under the alternative (complex) model.
#'
#' Under the null hypothesis and regularity conditions, \eqn{\Lambda}
#' asymptotically follows a chi-squared distribution with degrees of freedom
#' equal to the difference in number of free parameters.
#'
#' @examples
#' \dontrun{
#' # Generate data
#' set.seed(42)
#' x <- rnorm(100, mean = 5, sd = 2)
#'
#' # Fit full model (estimate both mu and sigma)
#' full <- fit(
#'   function(mu, log_sigma) {
#'     sigma <- exp(log_sigma)
#'     loglik_normal(mu, sigma, x)
#'   },
#'   params = c(mu = 0, log_sigma = 0)
#' )
#'
#' # Fit null model (mu fixed at 0)
#' null <- fit(
#'   function(log_sigma) {
#'     sigma <- exp(log_sigma)
#'     loglik_normal(0, sigma, x)  # mu = 0
#'   },
#'   params = c(log_sigma = 0)
#' )
#'
#' # Likelihood ratio test
#' test <- lrt(null, full)
#' test
#' pval(test)
#' is_significant_at(test, 0.05)
#'
#' # Can also use raw log-likelihoods
#' lrt(null_loglik = -150, alt_loglik = -140, df = 1)
#' }
#'
#' @importFrom stats pchisq
#' @export
lrt <- function(null_model, alt_model, df = NULL) {
  # Extract log-likelihoods
  if (inherits(null_model, "femtofit")) {
    null_ll <- null_model$loglik
    null_np <- length(coef(null_model))
  } else if (is.numeric(null_model)) {
    null_ll <- null_model
    null_np <- NULL
  } else {
    stop("null_model must be a femtofit object or numeric log-likelihood")
  }

  if (inherits(alt_model, "femtofit")) {
    alt_ll <- alt_model$loglik
    alt_np <- length(coef(alt_model))
  } else if (is.numeric(alt_model)) {
    alt_ll <- alt_model
    alt_np <- NULL
  } else {
    stop("alt_model must be a femtofit object or numeric log-likelihood")
  }

  # Compute degrees of freedom if not provided
  if (is.null(df)) {
    if (!is.null(null_np) && !is.null(alt_np)) {
      df <- alt_np - null_np
      if (df <= 0) {
        stop("Alternative model must have more parameters than null model")
      }
    } else {
      stop("df must be specified when using raw log-likelihoods")
    }
  }

  # Compute LRT statistic
  stat <- -2 * (null_ll - alt_ll)
  p.value <- pchisq(stat, df = df, lower.tail = FALSE)

  structure(
    list(
      stat = stat,
      p.value = p.value,
      dof = df,
      null_loglik = null_ll,
      alt_loglik = alt_ll
    ),
    class = c("likelihood_ratio_test", "hypothesis_test")
  )
}


#' Wald Test for Model Parameters
#'
#' Performs Wald tests for individual parameters or linear combinations
#' of parameters in a fitted model.
#'
#' @param object A `femtofit` object, or for the generic: the parameter estimate.
#' @param parm Parameter name(s) or indices to test. If NULL, tests all parameters.
#' @param null_value The null hypothesis value(s). Default is 0.
#' @param ... Additional arguments.
#'
#' @return For single parameter tests, a `wald_test` object (also inherits from
#'   `hypothesis_test`) containing:
#'   \describe{
#'     \item{stat}{The Wald chi-squared statistic}
#'     \item{p.value}{P-value from chi-squared(1) distribution}
#'     \item{dof}{Degrees of freedom (1 for single parameter)}
#'     \item{z}{The z-score}
#'     \item{estimate}{Parameter estimate}
#'     \item{se}{Standard error}
#'     \item{null_value}{The null hypothesis value}
#'   }
#'
#'   For multiple parameters, returns a list of wald_test objects.
#'
#' @details
#' The Wald test statistic for a single parameter is:
#' \deqn{W = \left(\frac{\hat{\theta} - \theta_0}{SE(\hat{\theta})}\right)^2}
#'
#' which follows a chi-squared distribution with 1 degree of freedom under H0.
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' x <- rnorm(100, mean = 5, sd = 2)
#'
#' result <- fit(
#'   function(mu, log_sigma) loglik_normal(mu, exp(log_sigma), x),
#'   params = c(mu = 0, log_sigma = 0)
#' )
#'
#' # Test if mu = 0
#' test <- wald_test(result, "mu")
#' pval(test)
#'
#' # Test if mu = 5
#' wald_test(result, "mu", null_value = 5)
#'
#' # Test all parameters
#' wald_test(result)
#' }
#'
#' @importFrom stats pchisq
#' @export
wald_test <- function(object, ...) {
  UseMethod("wald_test")
}


#' @describeIn wald_test Wald test for femtofit objects
#' @export
wald_test.femtofit <- function(object, parm = NULL, null_value = 0, ...) {
  cf <- coef(object)
  ses <- se(object)

  # Select parameters to test
  if (is.null(parm)) {
    parm <- names(cf)
  } else if (is.numeric(parm)) {
    parm <- names(cf)[parm]
  }

  # Recycle null_value if needed
  if (length(null_value) == 1) {
    null_value <- rep(null_value, length(parm))
  }
  names(null_value) <- parm

  # Compute Wald tests
  results <- lapply(parm, function(p) {
    est <- cf[p]
    se_val <- ses[p]
    null_val <- null_value[p]

    z <- (est - null_val) / se_val
    stat <- z^2
    p.value <- pchisq(stat, df = 1, lower.tail = FALSE)

    structure(
      list(
        stat = unname(stat),
        p.value = unname(p.value),
        dof = 1,
        z = unname(z),
        estimate = unname(est),
        se = unname(se_val),
        null_value = unname(null_val),
        parameter = p
      ),
      class = c("wald_test", "hypothesis_test")
    )
  })

  names(results) <- parm

  # Return single test or list
  if (length(results) == 1) {
    results[[1]]
  } else {
    results
  }
}


#' @describeIn wald_test Wald test from raw estimate and SE
#' @param se Standard error (for default method)
#' @export
wald_test.default <- function(object, se, null_value = 0, ...) {
  estimate <- object
  z <- (estimate - null_value) / se
  stat <- z^2
  p.value <- pchisq(stat, df = 1, lower.tail = FALSE)

  structure(
    list(
      stat = stat,
      p.value = p.value,
      dof = 1,
      z = z,
      estimate = estimate,
      se = se,
      null_value = null_value
    ),
    class = c("wald_test", "hypothesis_test")
  )
}


#' Extract p-value from hypothesis test
#'
#' @param x A hypothesis test object
#' @param ... Additional arguments (ignored)
#' @return The p-value
#' @export
pval <- function(x, ...) {
  UseMethod("pval")
}


#' @describeIn pval p-value for hypothesis tests
#' @export
pval.hypothesis_test <- function(x, ...) {
  x$p.value
}


#' Extract test statistic from hypothesis test
#'
#' @param x A hypothesis test object
#' @param ... Additional arguments (ignored)
#' @return The test statistic
#' @export
test_stat <- function(x, ...) {
  UseMethod("test_stat")
}


#' @describeIn test_stat Test statistic for hypothesis tests
#' @export
test_stat.hypothesis_test <- function(x, ...) {
  x$stat
}


#' Extract degrees of freedom from hypothesis test
#'
#' @param x A hypothesis test object
#' @param ... Additional arguments (ignored)
#' @return Degrees of freedom
#' @export
dof <- function(x, ...) {
  UseMethod("dof")
}


#' @describeIn dof Degrees of freedom for hypothesis tests
#' @export
dof.hypothesis_test <- function(x, ...) {
  x$dof
}


#' Check if test is significant at given level
#'
#' @param x A hypothesis test object
#' @param alpha Significance level (default 0.05)
#' @param ... Additional arguments (ignored)
#' @return Logical indicating significance
#' @export
is_significant_at <- function(x, alpha = 0.05, ...) {
  UseMethod("is_significant_at")
}


#' @describeIn is_significant_at Significance check for hypothesis tests
#' @export
is_significant_at.hypothesis_test <- function(x, alpha = 0.05, ...) {
  pval(x) < alpha
}


#' Print method for likelihood ratio test
#'
#' @param x A likelihood_ratio_test object
#' @param ... Additional arguments (ignored)
#' @export
print.likelihood_ratio_test <- function(x, ...) {
  cat("Likelihood Ratio Test\n")
  cat("---------------------\n")
  cat("LRT statistic:", format(x$stat, digits = 4), "\n")
  cat("Degrees of freedom:", x$dof, "\n")
  cat("P-value:", format(x$p.value, digits = 4), "\n")
  cat("\n")
  cat("Null log-likelihood:", format(x$null_loglik, digits = 4), "\n")
  cat("Alt log-likelihood:", format(x$alt_loglik, digits = 4), "\n")

  sig_level <- if (x$p.value < 0.001) "***" else if (x$p.value < 0.01) "**"
               else if (x$p.value < 0.05) "*" else ""
  if (nchar(sig_level) > 0) {
    cat("\nSignificance:", sig_level, "\n")
  }

  invisible(x)
}


#' Print method for Wald test
#'
#' @param x A wald_test object
#' @param ... Additional arguments (ignored)
#' @export
print.wald_test <- function(x, ...) {
  cat("Wald Test\n")
  cat("---------\n")

  if (!is.null(x$parameter)) {
    cat("Parameter:", x$parameter, "\n")
  }

  cat("Estimate:", format(x$estimate, digits = 4), "\n")
  cat("Std. Error:", format(x$se, digits = 4), "\n")
  cat("Null value:", x$null_value, "\n")
  cat("\n")
  cat("z-score:", format(x$z, digits = 4), "\n")
  cat("Wald statistic:", format(x$stat, digits = 4), "\n")
  cat("P-value:", format(x$p.value, digits = 4), "\n")

  sig_level <- if (x$p.value < 0.001) "***" else if (x$p.value < 0.01) "**"
               else if (x$p.value < 0.05) "*" else ""
  if (nchar(sig_level) > 0) {
    cat("\nSignificance:", sig_level, "\n")
  }

  invisible(x)
}
