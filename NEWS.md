`fwb` News and Updates
======

# `fwb` 0.5.1

* Fixed a bug where computing confidence intervals would yield an error about unused arguments for R versions prior to 4.5.0. Thanks to @vincentarelbundock. (#6)

# `fwb` 0.5.0

* Added a new confidence interval type for `confint()`, `fwb.ci()`, and `summary()`: `"wald"`, for Wald-type confidence intervals that don't correct for any bias.

* When p-values are requested in `summary()`, they are now based on inverting the confidence interval. This ensures hypothesis testing using the confidence interval and using p-values yield the same conclusion. Previously, they were based on inverting the Wald confidence interval only (i.e., a standard z-test).

* The null value of the estimates for the hypothesis tests in `summary()` can now be supplied using the `null` argument.

* Simultaneous inference via the sup-t confidence band and its inversion are now supported by `summary()` and `confint()` by setting `simultaneous = TRUE`. This is only supported for percentile and Wald confidence intervals (and the latter requires the `mvtnorm` package to be installed).

* Added new function `fwb.array()` to extract the bootstrap weights from an `fwb` object.

* Confidence intervals can be suppressed in `summary()` by setting `conf = 0`.

* Fixed a bug in `confint()`, `fwb.ci()`, and `summary()` where the confidence level could only be as low as .5. Now levels as low as just above 0 are allowed, except for when computing simultaneous Wald confidence intervals.

* BCa confidence intervals are computed faster in `confint()` and `summary()`. These functions no longer use `fwb.ci()` internally.

* Added a new `tidy()` method for `summary.fwb` objects.

# `fwb` 0.4.0

* Added a suite of new functions for computing weighted statistic and transformations that automatically incorporate the bootstrap weights. These include `w_mean()`, `w_var()`, `w_sd()`, `w_quantile()`, and `w_median()` for computing weighted means, variances, standard deviations, quantiles, and medians; `w_cov()` and `w_cor()` for computing weighted covariance and correlation matrices, and `w_std()`, `w_scale()`, and `w_center()` for transforming variables by standardizing, scaling, and centering using weighted statistics. These work when called inside the function supplied to the `statistic` argument of `fwb()` or inside the model that is supplied to `vcovFWB()`.

* Improved some error messages.

* Fixed a bug in `print.fwbci()` due to incorrect ordering of the intervals which led them to be printed with incorrect labels. These have been corrected and printing is a little prettier. Thanks to Katya Zelevinsky.

* Added `coef()` and `vcov()` methods for `fwb` objects.

* Documentation and vignette updates.

* Added new tests.

# `fwb` 0.3.0

* Added a new `confint()` method for `fwb` objects.

* Added a new `strata` argument to `fwb()` to perform stratified bootstrapping within levels of a stratification variable.

* Added a new `drop0` argument to `fwb()` to drop all units with weights of 0 in each bootstrap iteration.

* Added a new `.coef` argument to `vcovFWB()`. A function can be supplied to extract a vector of coefficients from the fitted model in each bootstrap iteration if the default (`stats::coef()`) doesn't return a numeric vector (e.g., for `nnet::multinom()` models). An error message is now thrown if `.coef` doesn't return a numeric vector.

* Added support for using `future` backend for `fwb()` by supplying `cl = "future"`. Thanks to Katya Zelevinsky for the suggestion.

* Added a new vignette on reproducibility and parallelization, which can be accessed at `vignette("fwb-rep")`.

* For `fwb()`, `simple` has a new default that is `TRUE` in most cases and `FALSE` when `wtype` is `"multinom"`. This should not affect results but will reduce memory use for large datasets by avoiding computing all bootstrap weights simultaneously. Note that when there is randomness in the `statistic` supplied to `fwb()`, the argument to `simple` affects whether BCa confidence intervals can be computed. See the reproducibility vignette mentioned above for details.

* A warning is now thrown when using `fwb()` with `simple = TRUE` with non-`NULL` `cl` when the random number generator kind is not `"L'Ecuyer-CMRG"`. Under these circumstances, results may not replicate and the BCa confidence interval will be inaccurate. See the reproducibility vignette mentioned above for details.

* Fixed a bug where the names of quantities produced by `fwb()` when `statistic` returns an unnamed vector were incorrect.

* When BCa confidence intervals are requested, an error is thrown if the number of bootstrap replications is smaller than the sample size.

* Documentation updates.

# `fwb` 0.2.0

* `fwb()` and `vcovFWB()` now take an additional argument, `wtype`, which specifies how the weights are drawn. The default, `"exp"` is still to draw weights from an $\text{Exp}(1)$ distribution but other options, namely `"multinom"` for multinomial integer weights (which reproduce `boot::boot()` results exactly), `"poisson"` for Poisson integer weights, and `"mammen"` for second-order accurate Mammen weights as recommended by Lihua Lei [here](https://x.com/lihua_lei_stat/status/1641538993090351106). (#4)

* New functions `set_fwb_wtype()` and `get_fwb_wtype()` allow one to set global defaults for the `wtype` argument of `fwb()` and vcovFWB()`.

# `fwb` 0.1.2

* Small updates and bug fixes.

# `fwb` 0.1.1

* Fixed bugs related to the `index` argument of various functions, including bugs when the estimated quantity is not given a name.

* Some error messages may be clearer.

# `fwb` 0.1.0

* First version!
