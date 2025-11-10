# Changelog

## `fwb` (development version)

- `wtype` can now be set to `"beta"` to sample weights from a
  \\\text{Beta}(1/2,3/2)\\ distribution and `"power"` to sample from a
  \\\text{Beta}(\sqrt{2} - 1, 1)\\ distribution as described by [Owen
  (2025)](https://doi.org/10.48550/arXiv.2508.10083).

- `drop0` can now be set to `NA` in
  [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md); this sets
  all weights of 0 to `NA` instead of removing those observations from
  the dataset.

- [`vcovFWB()`](https://ngreifer.github.io/fwb/reference/vcovFWB.md) now
  accepts `drop0` to control how to treat units with weights of 0.

- In [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md),
  `simple` can now be set to `TRUE` with `wtype = "multinom"`.
  `simple = FALSE` is still the default with `wtype = "multinom"` to
  maintain comparability with
  [`boot::boot()`](https://rdrr.io/pkg/boot/man/boot.html).

- [`vcovFWB()`](https://ngreifer.github.io/fwb/reference/vcovFWB.md) has
  improved support for `coxph` objects from *survival*.

- Added code of conduct to README.

- New tests.

## `fwb` 0.5.1

CRAN release: 2025-09-19

- Fixed a bug where computing confidence intervals would yield an error
  about unused arguments for R versions prior to 4.5.0. Thanks to
  [@vincentarelbundock](https://github.com/vincentarelbundock).
  ([\#6](https://github.com/ngreifer/fwb/issues/6))

## `fwb` 0.5.0

CRAN release: 2025-07-08

- Added a new confidence interval type for
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`fwb.ci()`](https://ngreifer.github.io/fwb/reference/fwb.ci.md), and
  [`summary()`](https://rdrr.io/r/base/summary.html): `"wald"`, for
  Wald-type confidence intervals that don’t correct for any bias.

- When p-values are requested in
  [`summary()`](https://rdrr.io/r/base/summary.html), they are now based
  on inverting the confidence interval. This ensures hypothesis testing
  using the confidence interval and using p-values yield the same
  conclusion. Previously, they were based on inverting the Wald
  confidence interval only (i.e., a standard z-test).

- The null value of the estimates for the hypothesis tests in
  [`summary()`](https://rdrr.io/r/base/summary.html) can now be supplied
  using the `null` argument.

- Simultaneous inference via the sup-t confidence band and its inversion
  are now supported by
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`confint()`](https://rdrr.io/r/stats/confint.html) by setting
  `simultaneous = TRUE`. This is only supported for percentile and Wald
  confidence intervals (and the latter requires the `mvtnorm` package to
  be installed).

- Added new function
  [`fwb.array()`](https://ngreifer.github.io/fwb/reference/fwb.array.md)
  to extract the bootstrap weights from an `fwb` object.

- Confidence intervals can be suppressed in
  [`summary()`](https://rdrr.io/r/base/summary.html) by setting
  `conf = 0`.

- Fixed a bug in [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`fwb.ci()`](https://ngreifer.github.io/fwb/reference/fwb.ci.md), and
  [`summary()`](https://rdrr.io/r/base/summary.html) where the
  confidence level could only be as low as .5. Now levels as low as just
  above 0 are allowed, except for when computing simultaneous Wald
  confidence intervals.

- BCa confidence intervals are computed faster in
  [`confint()`](https://rdrr.io/r/stats/confint.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html). These functions no
  longer use
  [`fwb.ci()`](https://ngreifer.github.io/fwb/reference/fwb.ci.md)
  internally.

- Added a new `tidy()` method for `summary.fwb` objects.

## `fwb` 0.4.0

CRAN release: 2025-06-11

- Added a suite of new functions for computing weighted statistic and
  transformations that automatically incorporate the bootstrap weights.
  These include
  [`w_mean()`](https://ngreifer.github.io/fwb/reference/w_mean.md),
  [`w_var()`](https://ngreifer.github.io/fwb/reference/w_mean.md),
  [`w_sd()`](https://ngreifer.github.io/fwb/reference/w_mean.md),
  [`w_quantile()`](https://ngreifer.github.io/fwb/reference/w_mean.md),
  and [`w_median()`](https://ngreifer.github.io/fwb/reference/w_mean.md)
  for computing weighted means, variances, standard deviations,
  quantiles, and medians;
  [`w_cov()`](https://ngreifer.github.io/fwb/reference/w_mean.md) and
  [`w_cor()`](https://ngreifer.github.io/fwb/reference/w_mean.md) for
  computing weighted covariance and correlation matrices, and
  [`w_std()`](https://ngreifer.github.io/fwb/reference/w_mean.md),
  [`w_scale()`](https://ngreifer.github.io/fwb/reference/w_mean.md), and
  [`w_center()`](https://ngreifer.github.io/fwb/reference/w_mean.md) for
  transforming variables by standardizing, scaling, and centering using
  weighted statistics. These work when called inside the function
  supplied to the `statistic` argument of
  [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) or inside
  the model that is supplied to
  [`vcovFWB()`](https://ngreifer.github.io/fwb/reference/vcovFWB.md).

- Improved some error messages.

- Fixed a bug in
  [`print.fwbci()`](https://ngreifer.github.io/fwb/reference/fwb.ci.md)
  due to incorrect ordering of the intervals which led them to be
  printed with incorrect labels. These have been corrected and printing
  is a little prettier. Thanks to Katya Zelevinsky.

- Added [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) methods for `fwb`
  objects.

- Documentation and vignette updates.

- Added new tests.

## `fwb` 0.3.0

CRAN release: 2025-03-03

- Added a new [`confint()`](https://rdrr.io/r/stats/confint.html) method
  for `fwb` objects.

- Added a new `strata` argument to
  [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) to perform
  stratified bootstrapping within levels of a stratification variable.

- Added a new `drop0` argument to
  [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) to drop all
  units with weights of 0 in each bootstrap iteration.

- Added a new `.coef` argument to
  [`vcovFWB()`](https://ngreifer.github.io/fwb/reference/vcovFWB.md). A
  function can be supplied to extract a vector of coefficients from the
  fitted model in each bootstrap iteration if the default
  ([`stats::coef()`](https://rdrr.io/r/stats/coef.html)) doesn’t return
  a numeric vector (e.g., for
  [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html)
  models). An error message is now thrown if `.coef` doesn’t return a
  numeric vector.

- Added support for using `future` backend for
  [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) by
  supplying `cl = "future"`. Thanks to Katya Zelevinsky for the
  suggestion.

- Added a new vignette on reproducibility and parallelization, which can
  be accessed at
  [`vignette("fwb-rep")`](https://ngreifer.github.io/fwb/articles/fwb-rep.md).

- For [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md),
  `simple` has a new default that is `TRUE` in most cases and `FALSE`
  when `wtype` is `"multinom"`. This should not affect results but will
  reduce memory use for large datasets by avoiding computing all
  bootstrap weights simultaneously. Note that when there is randomness
  in the `statistic` supplied to
  [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md), the
  argument to `simple` affects whether BCa confidence intervals can be
  computed. See the reproducibility vignette mentioned above for
  details.

- A warning is now thrown when using
  [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) with
  `simple = TRUE` with non-`NULL` `cl` when the random number generator
  kind is not `"L'Ecuyer-CMRG"`. Under these circumstances, results may
  not replicate and the BCa confidence interval will be inaccurate. See
  the reproducibility vignette mentioned above for details.

- Fixed a bug where the names of quantities produced by
  [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) when
  `statistic` returns an unnamed vector were incorrect.

- When BCa confidence intervals are requested, an error is thrown if the
  number of bootstrap replications is smaller than the sample size.

- Documentation updates.

## `fwb` 0.2.0

CRAN release: 2023-12-07

- [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) and
  [`vcovFWB()`](https://ngreifer.github.io/fwb/reference/vcovFWB.md) now
  take an additional argument, `wtype`, which specifies how the weights
  are drawn. The default, `"exp"` is still to draw weights from an
  \\\text{Exp}(1)\\ distribution but other options, namely `"multinom"`
  for multinomial integer weights (which reproduce
  [`boot::boot()`](https://rdrr.io/pkg/boot/man/boot.html) results
  exactly), `"poisson"` for Poisson integer weights, and `"mammen"` for
  second-order accurate Mammen weights as recommended by Lihua Lei
  [here](https://x.com/lihua_lei_stat/status/1641538993090351106).
  ([\#4](https://github.com/ngreifer/fwb/issues/4))

- New functions
  [`set_fwb_wtype()`](https://ngreifer.github.io/fwb/reference/set_fwb_wtype.md)
  and
  [`get_fwb_wtype()`](https://ngreifer.github.io/fwb/reference/set_fwb_wtype.md)
  allow one to set global defaults for the `wtype` argument of
  [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) and
  vcovFWB()\`.

## `fwb` 0.1.2

CRAN release: 2023-10-02

- Small updates and bug fixes.

## `fwb` 0.1.1

CRAN release: 2022-10-26

- Fixed bugs related to the `index` argument of various functions,
  including bugs when the estimated quantity is not given a name.

- Some error messages may be clearer.

## `fwb` 0.1.0

CRAN release: 2022-09-19

- First version!
