*fwb* News and Updates
======

# *fwb* 0.7.0

### New Features

* Added a new confidence interval type for `confint()`, `fwb.ci()`, and `summary()`: `"cheap"`, for the "cheap" confidence interval described by [Lam (2022)](https://arxiv.org/abs/2202.00090). This interval maintains nominal coverage even with few bootstrap replications (as low as one!).

### Breaking changes

* **Bootstrap weights now depend only on the seed.** Previously, when `simple = TRUE`, the weights were drawn inside each bootstrap replication and so depended on how the work was distributed: results changed with `cl`, with the number of workers, and with `verbose`, and reproducing a parallel run required `set.seed(###, "L'Ecuyer-CMRG")` (or, with a `cluster` object, `parallel::clusterSetRNGStream()` instead of `set.seed()`). `fwb()` and `vcovFWB()` now reserve one random number stream per replicate before any work is dispatched, and each replication draws its weights from its own stream. A plain call to `set.seed()` is all that is needed, whatever `cl`, `verbose`, or `RNGkind()` are set to, and `parallel::clusterSetRNGStream()` no longer has any effect. See `vignette("fwb-rep")`.

* **This changes the bootstrap estimates when `simple = TRUE`** (the default for every `wtype` except `"multinom"`). Analyses that need to reproduce results from *fwb* 0.6.0 or earlier should install that version (e.g., using `pak::pak("fwb@0.6.0")`). Results with `simple = FALSE` are unchanged, so `wtype = "multinom"` continues to match `boot::boot(., stype = "f")` exactly. Note also that `simple = TRUE` and `simple = FALSE` no longer give the same estimates as each other, even for the weight types where they previously did.

* `<fwb>` objects no longer store the `cl` argument, since nothing about the backend is needed to recover the weights any more. Previously, a `cluster` object kept in the result made the object unusable after `parallel::stopCluster()` or after being saved and reloaded.

* `vignette("fwb-rep")` has been rewritten around the new behavior. It no longer describes the several cases that used to require different treatment, and it explains how to reproduce results from earlier versions.

* **The multinomial bootstrap is now exhaustive when it can be.** There are only $n^n$ distinct multinomial bootstrap samples of $n$ units, so when `R` is at least that many, `fwb()` and `vcovFWB()` now use each of them exactly once instead of sampling, and set `R` to that number (reporting this in a message). The estimates are then the exact bootstrap distribution rather than a sample from it, and they no longer depend on the seed. This applies with `strata` (where the count is $\prod_s n_s^{n_s}$) and with `cluster` (where it is computed over clusters). It is the one respect in which `wtype = "multinom"` no longer matches `boot::boot(., stype = "f")`, which always samples; it reaches only very small samples, as $n = 6$ already has 46656 distinct samples.

### Other fixes

* Simultaneous inference with `ci.type = "perc"` is much faster and now exact. Previously the confidence level was found by numerical search, which could return a band whose coverage fell short of the one requested; it is now computed directly.

* Simultaneous inference now produces an error rather than `NA` when some of the bootstrap estimates are `NA` or non-finite.

* BCa confidence intervals and `fwb.array()` now work in every case. Previously, when `simple = TRUE` and the function supplied to `statistic` involved a random element, the weights could not be recovered: `fwb.array()` warned and returned weights that were not the ones used, and `fwb.ci()` and `summary()` refused to compute BCa intervals. Recovering the weights no longer involves replaying the random number stream, so these restrictions are gone.

* Fixed a bug in `fwb.array()` where the bootstrap weights were computed as if `cluster` had not been supplied, returning weights unrelated to those used in the cluster bootstrap.

* Fixed a bug in `fwb.array()` where `strata` was ignored unless it had been supplied as a `factor`, which also gave incorrect BCa confidence intervals for stratified bootstraps.

* Fixed a bug in `fwb.array()` where the bootstrap weights were incorrect when `wtype = "multinom"` and `simple = TRUE`.

* Fixed a bug in `fwb()` where stratified bootstrapping failed with an error when `simple = TRUE` (its default for every `wtype` except `"multinom"`).

* Fixed a bug in `vcovFWB()` that could give incorrect standard errors for `lm` models with a weight type that allows zeroes.

* Weight generation with `wtype = "multinom"` is faster, by up to 2x for small samples. Estimates are unchanged, so results still match `boot::boot(., stype = "f")` exactly.

* `fwb()` and `vcovFWB()` now throw an error when there is only one unit, or only one cluster, to resample. Previously `fwb()` failed with an obscure message for some weight types and returned all-zero standard errors for others, and `vcovFWB()` returned a matrix of zeros.

* `fwb()` now checks that the function supplied to `statistic` can accept the dataset and the weights, and that its arguments do not share a name with one of `fwb()`'s own (which would prevent them from ever being supplied through `...`). Functions with a `...` argument are exempt, since they can accept whatever they are given; this keeps packages that wrap `statistic`, such as [*progressify*](https://CRAN.R-project.org/package=progressify), working unchanged.

* Confidence interval types that cannot be computed from the available number of bootstrap replications now say so. Previously, with `R` below what a type needs, `"wald"` and `"norm"` returned `NA` limits without comment and the other types produced errors from deep inside the calculation ("subscript out of bounds", "missing value where TRUE/FALSE needed"). Only `"cheap"` can be computed from a single replication; the others require at least two.

* `fwb.ci()` now computes the interval types it can and skips the rest with a warning, rather than failing outright. Previously `type = "all"` at a small `R` failed on the first type it could not compute, returning nothing.

* `summary()` and `confint()` now enforce BCa's requirements, as `fwb.ci()` always has. Previously, requesting `ci.type = "bca"` with fewer bootstrap replications than units, or with clusters, produced an interval from an underdetermined calculation without any indication that it was not usable.

* Fixed a bug in `vcovFWB()` where an error was produced for any model whose fit dropped rows due to missing values in the model variables (e.g., under the default `na.action`). The weights are now aligned with the rows the model used before the model is re-fit.

* `vcovFWB()` is faster when `drop0` is `TRUE` or `NA` for `lm` and `glm` models. Previously, these settings silently disabled the fast re-fitting paths used for these models, so every replicate was re-fit with a full `update()` call.

* Fixed a bug in `vcovFWB()` where `drop0 = NA` produced an error ("cannot find valid starting values") for `glm` models. All three values of `drop0` now give identical results for `lm` and `glm` models, as they should.

* `vcovFWB()` no longer fails for models that reject weights of 0 (e.g., `survival::coxph()`) while it determines how to re-fit them.

* Fixed a bug in `vcovFWB()` where `fix = TRUE` produced an error whenever the covariance matrix it was meant to correct was not positive semi-definite.

* Fixed a bug in `vcovFWB()` where an error was produced in a session that had not yet used the random number generator (e.g., in a fresh `Rscript` process).

* Fixed a bug in `vcovFWB()` where `drop0 = NA` set the weights of units with nonzero weights to `NA` rather than those with weights of 0.

* Fixed a bug in `summary()` where simultaneous p-values were assigned to the wrong estimates when `ci.type = "wald"` and one of the estimates had a variance of 0.

* Fixed a bug in `w_scale()` (and `w_std()` with `center = FALSE`) where the variable was scaled by its uncentered second moment rather than by its weighted standard deviation as documented, when weights were supplied.

* Fixed an incorrect error when using `confint.fwb()` with incorrectly specified `parm`.

* `w_std()`, `w_scale()`, and `w_center()` now always return a vector as long as their input. Previously, units with weights of 0 or with missing values were dropped from the output, which made these functions fail inside a model formula (as demonstrated in `help("w_std")`) whenever any weight was 0, i.e., with `wtype = "multinom"` or `"poisson"`.

* The `w_*()` functions can now be used inside `statistic` when `cl` is a `cluster` object; previously this failed with an error about the function not being found, because *fwb* was not attached on the workers.

* `w_cor()` gains an `na.rm` argument, matching the rest of the `w_*()` functions.

* `vcovFWB()` now requires `R` to be greater than 1; previously `R = 0` produced an unrelated error from `stats::cov()`.

* `vcovFWB()` no longer forwards `...` to the model refitting function; as documented, it is ignored.

* Documentation fixes.

* New tests.

# *fwb* 0.6.0

* `wtype` can now be set to `"beta"` to sample weights from a $\text{Beta}(1/2,3/2)$ distribution or `"power"` to sample from a $\text{Beta}(\sqrt{2} - 1, 1)$ distribution as described by [Owen (2025)](https://doi.org/10.48550/arXiv.2508.10083).

* `drop0` can now be set to `NA` in `fwb()`; this sets all weights of 0 to `NA` instead of removing those observations from the dataset.

* `vcovFWB()` now accepts `drop0` to control how to treat units with weights of 0.

* In `fwb()`, `simple` can now be set to `TRUE` with `wtype = "multinom"`. `simple = FALSE` is still the default with `wtype = "multinom"` to maintain comparability with `boot::boot()`.

* `vcovFWB()` has improved support for `coxph` objects from *survival*.

* In `fwb()`, `verbose` is now `FALSE` by default when parallelization is used.

* [*arg*](https://ngreifer.github.io/arg/) is now used for errors and warning messages.

* Updated `vignette("fwb-rep")` to discuss reproducibility with the `verbose` argument.

* Added code of conduct to README.

* New tests.

* Documentation updates.

# *fwb* 0.5.1

* Fixed a bug where computing confidence intervals would yield an error about unused arguments for R versions prior to 4.5.0. Thanks to @vincentarelbundock. (#6)

# *fwb* 0.5.0

* Added a new confidence interval type for `confint()`, `fwb.ci()`, and `summary()`: `"wald"`, for Wald-type confidence intervals that don't correct for any bias.

* When p-values are requested in `summary()`, they are now based on inverting the confidence interval. This ensures hypothesis testing using the confidence interval and using p-values yield the same conclusion. Previously, they were based on inverting the Wald confidence interval only (i.e., a standard z-test).

* The null value of the estimates for the hypothesis tests in `summary()` can now be supplied using the `null` argument.

* Simultaneous inference via the sup-t confidence band and its inversion are now supported by `summary()` and `confint()` by setting `simultaneous = TRUE`. This is only supported for percentile and Wald confidence intervals (and the latter requires the `mvtnorm` package to be installed).

* Added new function `fwb.array()` to extract the bootstrap weights from an `<fwb>` object.

* Confidence intervals can be suppressed in `summary()` by setting `conf = 0`.

* Fixed a bug in `confint()`, `fwb.ci()`, and `summary()` where the confidence level could only be as low as .5. Now levels as low as just above 0 are allowed, except for when computing simultaneous Wald confidence intervals.

* BCa confidence intervals are computed faster in `confint()` and `summary()`. These functions no longer use `fwb.ci()` internally.

* Added a new `tidy()` method for `summary.fwb` objects.

# *fwb* 0.4.0

* Added a suite of new functions for computing weighted statistic and transformations that automatically incorporate the bootstrap weights. These include `w_mean()`, `w_var()`, `w_sd()`, `w_quantile()`, and `w_median()` for computing weighted means, variances, standard deviations, quantiles, and medians; `w_cov()` and `w_cor()` for computing weighted covariance and correlation matrices, and `w_std()`, `w_scale()`, and `w_center()` for transforming variables by standardizing, scaling, and centering using weighted statistics. These work when called inside the function supplied to the `statistic` argument of `fwb()` or inside the model that is supplied to `vcovFWB()`.

* Improved some error messages.

* Fixed a bug in `print.fwbci()` due to incorrect ordering of the intervals which led them to be printed with incorrect labels. These have been corrected and printing is a little prettier. Thanks to Katya Zelevinsky.

* Added `coef()` and `vcov()` methods for `<fwb>` objects.

* Documentation and vignette updates.

* Added new tests.

# *fwb* 0.3.0

* Added a new `confint()` method for `<fwb>` objects.

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

# *fwb* 0.2.0

* `fwb()` and `vcovFWB()` now take an additional argument, `wtype`, which specifies how the weights are drawn. The default, `"exp"` is still to draw weights from an $\text{Exp}(1)$ distribution but other options, namely `"multinom"` for multinomial integer weights (which reproduce `boot::boot()` results exactly), `"poisson"` for Poisson integer weights, and `"mammen"` for second-order accurate Mammen weights as recommended by Lihua Lei [here](https://x.com/lihua_lei_stat/status/1641538993090351106). (#4)

* New functions `set_fwb_wtype()` and `get_fwb_wtype()` allow one to set global defaults for the `wtype` argument of `fwb()` and vcovFWB()`.

# *fwb* 0.1.2

* Small updates and bug fixes.

# *fwb* 0.1.1

* Fixed bugs related to the `index` argument of various functions, including bugs when the estimated quantity is not given a name.

* Some error messages may be clearer.

# *fwb* 0.1.0

* First version!
