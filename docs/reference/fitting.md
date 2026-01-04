# Statistical model fitting with automatic differentiation

This module provides the
[`fit()`](https://queelius.github.io/femtograd/reference/fit.md)
function for maximum likelihood estimation and the `femtofit` class for
representing fitted models. The interface follows standard R
conventions, implementing base R generics like
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html), etc.
