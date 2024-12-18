`cobalt` News and Updates
======

# `fwb` (development version)

* Add a new `.coef` argument to `vcovFWB()`. A function can be supplied to extract a vector of coefficients from the fitted model in each bootstrap iteration if the default (`stats::coef()`) doesn't return a numeric vector (e.g., for `nnet::multinom()` models). An error message is now thrown if `.coef` doesn't return a numeric vector.

* Added support for using `future` backend for `fwb()` by supplying `cl = "future"`. Thanks to Katya Zelevinsky for the suggestion.

* A warning is now thrown when using `fwb()` with `simple = TRUE` with non-`NULL` `cl` when the random number generator kind is not `"L'Ecuyer-CMRG"`. Under these circumstances, results may not replicate and the BCa confidence interval will be inaccurate.

* Fixed a bug where the names of quantities produced by `fwb()` when `statistic` returns an unnamed vector were incorrect.

# `fwb` 0.2.0

* `fwb()` and `vcovFWB()` now take an additional argument, `wtype`, which specifies how the weights are drawn. The default, `"exp"` is still to draw weights from an $\text{Exp}(1)$ distribution but other options, namely `"multinom"` for multinomial integer weights (which reproduce `boot::boot()` results exactly), `"poisson"` for Poisson integer weights, and `"mammen"` for second-order accurate Mammen weights as recommended by Lihua Lei [here](https://twitter.com/lihua_lei_stat/status/1641538993090351106). (#4)

* New functions `set_fwb_wtype()` and `get_fwb_wtype()` allow one to set global defaults for the `wtype` argument of `fwb()` and vcovFWB()`.

# `fwb` 0.1.2

* Small updates and bug fixes.

# `fwb` 0.1.1

* Fixed bugs related to the `index` argument of various functions, including bugs when the estimated quantity is not given a name.

* Some error messages may be clearer.

# `fwb` 0.1.0

* First version!
