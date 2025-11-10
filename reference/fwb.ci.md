# Fractional Weighted Bootstrap Confidence Intervals

`fwb.ci()` generates several types of equi-tailed two-sided
nonparametric confidence intervals. These include the normal
approximation, the basic bootstrap interval, the percentile bootstrap
interval, the bias-corrected percentile bootstrap interval, and the
bias-correct and accelerated (BCa) bootstrap interval.

## Usage

``` r
fwb.ci(
  fwb.out,
  conf = 0.95,
  type = "bc",
  index = 1L,
  h = base::identity,
  hinv = base::identity,
  ...
)

# S3 method for class 'fwbci'
print(x, hinv = NULL, ...)
```

## Arguments

- fwb.out:

  an `fwb` object; the output of a call to
  [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md).

- conf:

  the desired confidence level. Default is .95 for 95% confidence
  intervals.

- type:

  the type of confidence interval desired. Allowable options include
  `"wald"` (Wald interval), `"norm"` (normal approximation), `"basic"`
  (basic interval), `"perc"` (percentile interval), `"bc"` (bias-correct
  percentile interval), and `"bca"` (BCa interval). More than one is
  allowed. Can also be `"all"` to request all of them. BCa intervals
  require that the number of bootstrap replications is larger than the
  sample size.

- index:

  the index of the position of the quantity of interest in `fwb.out$t0`
  if more than one was specified in
  [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md). Only one
  value is allowed at a time. By default the first statistic is used.

- h:

  a function defining a transformation. The intervals are calculated on
  the scale of `h(t)` and the inverse function `hinv` applied to the
  resulting intervals. It must be a function of one variable only and
  for a vector argument, it must return a vector of the same length.
  Default is the identity function.

- hinv:

  a function, like `h`, which returns the inverse of `h`. It is used to
  transform the intervals calculated on the scale of `h(t)` back to the
  original scale. The default is the identity function. If `h` is
  supplied but `hinv` is not, then the intervals returned will be on the
  transformed scale.

- ...:

  ignored

- x:

  an `fwbci` object; the output of a call to `fwb.ci()`.

## Value

An `fwbci` object, which inherits from `bootci` and has the following
components:

- R:

  the number of bootstrap replications in the original call to
  [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md).

- t0:

  the observed value of the statistic on the same scale as the intervals
  (i.e., after applying `h` and then `hinv`.

- call:

  the call to `fwb.ci()`.

There will be additional components named after each confidence interval
type requested. For `"wald"` and `"norm"`, this is a matrix with one row
containing the confidence level and the two confidence interval limits.
For the others, this is a matrix with one row containing the confidence
level, the indices of the two order statistics used in the calculations,
and the confidence interval limits.

## Details

`fwb.ci()` functions similarly to
[`boot::boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html) in that
it takes in a bootstrapped object and computes confidence intervals.
This interface is a bit old-fashioned, but was designed to mimic that of
`boot.ci()`. For a more modern interface, see
[`summary.fwb()`](https://ngreifer.github.io/fwb/reference/summary.fwb.md).

The bootstrap intervals are defined as follows, with \\\alpha =\\ 1 -
`conf`, \\t_0\\ the estimate in the original sample, \\\hat{t}\\ the
average of the bootstrap estimates, \\s_t\\ the standard deviation of
the bootstrap estimates, \\t^{(i)}\\ the set of ordered estimates with
\\i\\ corresponding to their quantile, and \\z\_\frac{\alpha}{2}\\ and
\\z\_{1-\frac{\alpha}{2}}\\ the upper and lower critical \\z\\ scores.

- `"wald"`:

  \$\$\left\[t_0 + s_t z\_\frac{\alpha}{2}, t_0 + s_t
  z\_{1-\frac{\alpha}{2}}\right\]\$\$ This method is valid when the
  statistic is normally distributed around the estimate.

- `"norm"` (normal approximation):

  \$\$\left\[2t_0 - \hat{t} + s_t z\_\frac{\alpha}{2}, 2t_0 - \hat{t} +
  s_t z\_{1-\frac{\alpha}{2}}\right\]\$\$

  This involves subtracting the "bias" (\\\hat{t} - t_0\\) from the
  estimate \\t_0\\ and using a standard Wald-type confidence interval.
  This method is valid when the statistic is normally distributed.

- `"basic"`:

  \$\$\left\[2t_0 - t^{(1-\frac{\alpha}{2})}, 2t_0 -
  t^{(\frac{\alpha}{2})}\right\]\$\$

- `"perc"` (percentile confidence interval):

  \$\$\left\[t^{(\frac{\alpha}{2})},
  t^{(1-\frac{\alpha}{2})}\right\]\$\$

- `"bc"` (bias-corrected percentile confidence interval):

  \$\$\left\[t^{(l)}, t^{(u)}\right\]\$\$

  \\l = \Phi\left(2z_0 + z\_\frac{\alpha}{2}\right)\\, \\u =
  \Phi\left(2z_0 + z\_{1-\frac{\alpha}{2}}\right)\\, where \\\Phi(.)\\
  is the normal cumulative density function (i.e.,
  [`pnorm()`](https://rdrr.io/r/stats/Normal.html)) and \\z_0 =
  \Phi^{-1}(q)\\ where \\q\\ is the proportion of bootstrap estimates
  less than the original estimate \\t_0\\. This is similar to the
  percentile confidence interval but changes the specific quantiles of
  the bootstrap estimates to use, correcting for bias in the original
  estimate. It is described in Xu et al. (2020). When \\t^0\\ is the
  median of the bootstrap distribution, the `"perc"` and `"bc"`
  intervals coincide.

- `"bca"` (bias-corrected and accelerated confidence interval):

  \$\$\left\[t^{(l)}, t^{(u)}\right\]\$\$

  \\l = \Phi\left(z_0 + \frac{z_0 +
  z\_\frac{\alpha}{2}}{1-a(z_0+z\_\frac{\alpha}{2})}\right)\\, \\u =
  \Phi\left(z_0 + \frac{z_0 +
  z\_{1-\frac{\alpha}{2}}}{1-a(z_0+z\_{1-\frac{\alpha}{2}})}\right)\\,
  using the same definitions as above, but with the additional
  acceleration parameter \\a\\, where \\a =
  \frac{1}{6}\frac{\sum{L^3}}{(\sum{L^2})^{3/2}}\\. \\L\\ is the
  empirical influence value of each unit, which is computed using the
  regression method described in
  [`boot::empinf()`](https://rdrr.io/pkg/boot/man/empinf.html). When
  \\a=0\\, the `"bca"` and `"bc"` intervals coincide. The acceleration
  parameter corrects for bias and skewness in the statistic. It can only
  be used when clusters are absent and the number of bootstrap
  replications is larger than the sample size. Note that BCa intervals
  cannot be requested when `simple = TRUE` and there is randomness in
  the `statistic` supplied to
  [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md).

Interpolation on the normal quantile scale is used when a non-integer
order statistic is required, as in
[`boot::boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html). Note
that unlike with
[`boot::boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html),
studentized confidence intervals (`type = "stud"`) are not allowed.

## Functions

- `print(fwbci)`: Print a bootstrap confidence interval

## See also

- [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) for
  performing the fractional weighted bootstrap

- [`get_ci()`](https://ngreifer.github.io/fwb/reference/get_ci.md) for
  extracting confidence intervals from an `fwbci` object

- [`summary.fwb()`](https://ngreifer.github.io/fwb/reference/summary.fwb.md)
  for producing clean output from
  [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) that
  includes confidence intervals calculated by `fwb.ci()`

- [`boot::boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html) for
  computing confidence intervals from the traditional bootstrap

- [`vcovFWB()`](https://ngreifer.github.io/fwb/reference/vcovFWB.md) for
  computing parameter estimate covariance matrices using the fractional
  weighted bootstrap

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

# Bias corrected percentile interval
bcci <- fwb.ci(fwb_out, index = "spontaneous",
               type = "bc")
bcci
#> BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
#>     Based on 199 bootstrap replicates     
#> 
#> CALL :
#> fwb.ci(fwb.out = fwb_out, type = "bc", index = "spontaneous")
#> 
#> Intervals :
#> 
#> Level BC Percentile 
#>  95%  (0.745, 1.638)
#> 
#> Calculations and Intervals on Original Scale
#> Some bias-corrected percentile intervals may be unstable

# Using `get_ci()` to extract confidence limits

get_ci(bcci)
#> $bc
#>         L         U 
#> 0.7453669 1.6378339 
#> 
#> attr(,"conf")
#> [1] 0.95

# Interval calculated on original (log odds) scale,
# then transformed by exponentiation to be on OR
fwb.ci(fwb_out, index = "induced", type = "norm",
       hinv = exp)
#> BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
#>     Based on 199 bootstrap replicates     
#> 
#> CALL :
#> fwb.ci(fwb.out = fwb_out, type = "norm", index = "induced", hinv = exp)
#> 
#> Intervals :
#> 
#> Level     Normal    
#>  95%  (1.029, 2.232)
#> 
#> Calculations on Original Scale but Intervals Transformed
```
