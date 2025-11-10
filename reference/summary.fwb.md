# Summarize `fwb` Output

[`summary()`](https://rdrr.io/r/base/summary.html) creates a regression
summary-like table that displays the bootstrap estimates, their
empirical standard errors, their confidence intervals, and, optionally,
p-values for tests against a null value.
[`confint()`](https://rdrr.io/r/stats/confint.html) produces just the
confidence intervals, and is called internally by
[`summary()`](https://rdrr.io/r/base/summary.html).

## Usage

``` r
# S3 method for class 'fwb'
summary(
  object,
  conf = 0.95,
  ci.type = "bc",
  p.value = FALSE,
  index = seq_len(ncol(object$t)),
  null = 0,
  simultaneous = FALSE,
  ...
)

# S3 method for class 'fwb'
confint(object, parm, level = 0.95, ci.type = "bc", simultaneous = FALSE, ...)
```

## Arguments

- object:

  an `fwb` object; the output of a call to
  [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md).

- conf, level:

  the desired confidence level. Default is .95 for 95% confidence
  intervals. Set to 0 to prevent calculation of confidence intervals.

- ci.type:

  the type of confidence interval desired. Allowable options include
  `"wald"` (Wald interval), `"norm"` (normal approximation), `"basic"`
  (basic interval), `"perc"` (percentile interval), `"bc"`
  (bias-corrected percentile interval), and `"bca"` (bias-corrected and
  accelerated \[BCa\] interval). Only one is allowed. BCa intervals
  require the number of bootstrap replications to be larger than the
  sample size. See
  [`fwb.ci()`](https://ngreifer.github.io/fwb/reference/fwb.ci.md) for
  details. The default is `"bc"`. Ignored if both `conf = 0` and
  `p.values = FALSE`.

- p.value:

  `logical`; whether to display p-values for the test that each
  parameter is equal to `null`. Default is `FALSE`. See Details.

- index, parm:

  the index or indices of the position of the quantity of interest if
  more than one was specified in
  [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md). Default is
  to display all quantities.

- null:

  `numeric`; when `p.value = TRUE`, the value of the estimate under the
  null hypothesis. Default is 0. Only one value can be supplied and it
  is applied to all tests.

- simultaneous:

  `logical`; whether to adjust confidence intervals and p-values to
  ensure correct familywise coverage/size across all specified
  estimates. See Details. Default is `FALSE` for standard pointwise
  intervals. `TRUE` is only allowed when `ci.type` is `"wald"` or
  `"perc"`.

- ...:

  ignored.

## Value

For [`summary()`](https://rdrr.io/r/base/summary.html), a `summary.fwb`
object, which is a matrix with the following columns:

- `Estimate`: the statistic estimated in the original sample

- `Std. Error`: the standard deviation of the bootstrap estimates

- `CI {L}%` and `CI {U}%`: the upper and lower confidence interval
  bounds computed using the argument to `ci.type` (only when `conf` is
  not 0).

- `z value`: when `p.value = TRUE` and `ci.type = "wald"`, the
  z-statistic for the test of the estimate against against `null`.

- `Pr(>|z|)`: when `p.value = TRUE`, the p-value for the test of the
  estimate against against `null`.

For [`confint()`](https://rdrr.io/r/stats/confint.html), a matrix with a
row for each statistic and a column for the upper and lower confidence
interval limits.

## Details

P-values are computed by inverting the confidence interval for each
parameter, i.e., finding the largest confidence level yielding a
confidence interval that excludes `null`, and taking the p-value to be
one minus that level. This ensures conclusions from tests based on the
p-value and whether the confidence interval contains the null value
always yield the same conclusion. Prior to version 0.5.0, all p-values
were based on inverting Wald confidence intervals, regardless of
`ci.type`.

Simultaneous confidence intervals are computed using the "sup-t"
confidence band, which involves modifying the confidence level so that
the intersection of all the adjusted confidence intervals contain the
whole parameter vector with the specified coverage. This will always be
less conservative than Bonferroni or Holm adjustment. See Olea and
Plagborg-Møller (2019) for details on implementation for Wald and
percentile intervals. Simultaneous p-values are computed by inverting
the simultaneous bands. Simultaneous inference is only allowed when
`ci.type` is `"wald"` or `"perc"` and `index` has length greater than 1.
When `ci.type = "wald"`, the mvtnorm package must be installed.

`tidy()` and [`print()`](https://rdrr.io/r/base/print.html) methods are
available for `summary.fwb` objects.

## References

Montiel Olea, J. L., & Plagborg-Møller, M. (2019). Simultaneous
confidence bands: Theory, implementation, and an application to SVARs.
*Journal of Applied Econometrics*, 34(1), 1–17.
[doi:10.1002/jae.2656](https://doi.org/10.1002/jae.2656)

## See also

[`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) for
performing the fractional weighted bootstrap;
[`fwb.ci()`](https://ngreifer.github.io/fwb/reference/fwb.ci.md) for
computing multiple confidence intervals for a single bootstrapped
quantity

## Examples

``` r
set.seed(123, "L'Ecuyer-CMRG")
data("infert")

fit_fun <- function(data, w) {
  fit <- glm(case ~ spontaneous + induced, data = data,
             family = "quasibinomial", weights = w)
  coef(fit)
}

fwb_out <- fwb(infert, fit_fun, R = 199,
               verbose = FALSE)

# Basic confidence interval for both estimates
summary(fwb_out, ci.type = "basic")
#>             Estimate Std. Error CI 2.5 % CI 97.5 %
#> (Intercept)  -1.7079     0.2647  -2.2118   -1.2074
#> spontaneous   1.1972     0.2185   0.7532    1.6483
#> induced       0.4181     0.1975   0.0331    0.8116

# Just for "induced" coefficient; p-values requested,
# no confidence intervals
summary(fwb_out, ci.type = "norm", conf = 0,
        index = "induced", p.value = TRUE)
#>         Estimate Std. Error Pr(>|z|)  
#> induced    0.418      0.197    0.035 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
