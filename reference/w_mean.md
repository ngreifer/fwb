# Calculate weighted statistics

These functions are designed to compute weighted statistics (mean,
variance, standard deviation, covariance, correlation, median and
quantiles) to perform weighted transformation (scaling, centering, and
standardization) using bootstrap weights. These automatically extract
bootstrap weights when called from within
[`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) to facilitate
computation of bootstrap statistics without needing to think to hard
about how to correctly incorporate weights.

## Usage

``` r
w_mean(x, w = NULL, na.rm = FALSE)

w_var(x, w = NULL, na.rm = FALSE)

w_sd(x, w = NULL, na.rm = FALSE)

w_cov(x, w = NULL, na.rm = FALSE)

w_cor(x, w = NULL)

w_quantile(
  x,
  w = NULL,
  probs = seq(0, 1, by = 0.25),
  na.rm = FALSE,
  names = TRUE,
  type = 7L,
  digits = 7L
)

w_median(x, w = NULL, na.rm = FALSE)

w_std(x, w = NULL, na.rm = TRUE, scale = TRUE, center = TRUE)

w_scale(x, w = NULL, na.rm = TRUE)

w_center(x, w = NULL, na.rm = TRUE)
```

## Arguments

- x:

  a numeric variable; for `w_cov()` and `w_cor()`, a numeric matrix.

- w:

  optional; a numeric vector of weights with length equal to the length
  of `x` (or number of rows if `x` is a matrix). If unspecified, will
  use bootstrapped weights when called from within
  [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) or
  [`vcovFWB()`](https://ngreifer.github.io/fwb/reference/vcovFWB.md) and
  unit weights (i.e., for unweighted estimates) otherwise.

- na.rm:

  logical; whether to exclude `NA` values in the weights and `x` when
  computing statistics. Default is `FALSE` for the weighted statistics
  (like with their unweighted counterparts) and `TRUE` for weighted
  transformations.

- probs, names, type, digits:

  see [`quantile()`](https://rdrr.io/r/stats/quantile.html). Only
  `type = 7` is allowed.

- scale, center:

  logical; whether to scale or center the variable.

## Value

`w_mean()`, `w_var()`, `w_sd()`, and `w_median()` return a numeric
scalar. `w_cov()` and `w_cor()` return numeric matrices. `w_quantile()`
returns a numeric vector of length equal to `probs`. `w_std()`,
`w_scale()`, and `w_center()` return numeric vectors of length equal to
the length of `x`.

## Details

These function automatically incorporate bootstrap weights when used
inside [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) or
[`vcovFWB()`](https://ngreifer.github.io/fwb/reference/vcovFWB.md). This
works because [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md)
and [`vcovFWB()`](https://ngreifer.github.io/fwb/reference/vcovFWB.md)
temporarily set an option with the environment of the function that
calls the estimation function with the sampled weights, and the `w_*()`
functions access the bootstrap weights from that environment, if any.
So, using, e.g., `w_mean()` inside the function supplied to the
`statistic` argument of
[`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md), computes the
weighted mean using the bootstrap weights. Using these functions outside
[`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) works just
like other functions that compute weighted statistics: when `w` is
supplied, the statistics are weighted, and otherwise they are
unweighted.

See below for how each statistic is computed.

### Weighted statistics

For all weighted statistics, the weights are first rescaled to sum to 1.
\\w\\ in the formulas below refers to these weights.

- `w_mean()`:

  Calculates the weighted mean as \$\$\bar{x}\_w = \sum_i{w_i x_i}\$\$
  This is the same as
  [`weighted.mean()`](https://rdrr.io/r/stats/weighted.mean.html).

- `w_var()`:

  Calculates the weighted variance as \$\$s^2\_{x,w} =
  \frac{\sum_i{w_i(x_i - \bar{x}\_w)^2}}{1 - \sum_i{w_i^2}}\$\$

- `w_sd()`:

  Calculates the weighted standard deviation as \$\$s\_{x,w} =
  \sqrt{s^2\_{x,w}}\$\$

- `w_cov()`:

  Calculates the weighted covariances as \$\$s\_{x,y,w} =
  \frac{\sum_i{w_i(x_i - \bar{x}\_w)(y_i - \bar{y}\_w)}}{1 -
  \sum_i{w_i^2}}\$\$ This is the same as
  [`cov.wt()`](https://rdrr.io/r/stats/cov.wt.html).

- `w_cor()`:

  Calculates the weighted correlation as \$\$r\_{x,y,w} =
  \frac{s\_{x,y,w}}{s\_{x,w}s\_{y,w}}\$\$ This is the same as
  [`cov.wt()`](https://rdrr.io/r/stats/cov.wt.html).

- `w_quantile()`:

  Calculates the weighted quantiles using linear interpolation of the
  weighted cumulative density function.

- `w_median()`:

  Calculates the weighted median as `w_quantile(., probs = .5)`.

### Weighted transformations

Weighted transformations use the weighted mean and/or standard deviation
to re-center or re-scale the given variable. In the formulas below,
\\x\\ is the original variable and \\w\\ are the weights.

- `w_center()`:

  Centers the variable at its (weighted) mean. \$\$x\_{c,w} = x -
  \bar{x}\_w\$\$

- `w_scale()`:

  Scales the variable by its (weighted) standard deviation. \$\$x\_{s,w}
  = x / s\_{x,w}\$\$

- `w_std()`:

  Centers and scales the variable by its (weighted) mean and standard
  deviation. \$\$x\_{z,w} = (x - \bar{x}\_w) / s\_{x,w}\$\$

`w_scale()` and `w_center()` are efficient wrappers to `w_std()` with
`center = FALSE` and `scale = FALSE`, respectively.

## See also

- [`mean()`](https://rdrr.io/r/base/mean.html) and
  [`weighted.mean()`](https://rdrr.io/r/stats/weighted.mean.html) for
  computing the unweighted and weighted mean

- [`var()`](https://rdrr.io/r/stats/cor.html) and
  [`sd()`](https://rdrr.io/r/stats/sd.html) for computing the unweighted
  variance and standard deviation

- [`median()`](https://rdrr.io/r/stats/median.html) and
  [`quantile()`](https://rdrr.io/r/stats/quantile.html) for computing
  the unweighted median and quantiles

- [`cov()`](https://rdrr.io/r/stats/cor.html),
  [`cor()`](https://rdrr.io/r/stats/cor.html), and
  [`cov.wt()`](https://rdrr.io/r/stats/cov.wt.html) for unweighted and
  weighted covariance and correlation matrices

- [`scale()`](https://rdrr.io/r/base/scale.html) for standardizing
  variables using arbitrary (by default, unweighted) centering and
  scaling factors

## Examples

``` r
# G-computation of average treatment effects using lalonde
# dataset

data("lalonde", package = "cobalt")

ate_est <- function(data, w) {
  fit <- lm(re78 ~ treat * (age + educ + married + race +
                              nodegree + re74 + re75),
            data = data, weights = w)

  p0 <- predict(fit, newdata = transform(data, treat = 0))
  p1 <- predict(fit, newdata = transform(data, treat = 1))

  # Weighted means using bootstrap weights
  m0 <- w_mean(p0)
  m1 <- w_mean(p1)

  c(m0 = m0, m1 = m1, ATE = m1 - m0)
}

set.seed(123, "L'Ecuyer-CMRG")
boot_est <- fwb(lalonde, statistic = ate_est,
                R = 199, verbose = FALSE)
summary(boot_est)
#>     Estimate Std. Error CI 2.5 % CI 97.5 %
#> m0      6296        344     5587      6978
#> m1      7371        974     5702      9688
#> ATE     1075        987     -672      3167

# Using `w_*()` data transformations inside a model
# supplied to vcovFWB():
fit <- lm(re78 ~ treat * w_center(age),
          data = lalonde)

lmtest::coeftest(fit, vcov = vcovFWB, R = 500)
#> 
#> t test of coefficients:
#> 
#>                     Estimate Std. Error t value Pr(>|t|)    
#> (Intercept)         6934.358    366.357 18.9279  < 2e-16 ***
#> treat               -436.121    682.525 -0.6390  0.52307    
#> w_center(age)         74.667     37.289  2.0024  0.04569 *  
#> treat:w_center(age)   21.710     72.651  0.2988  0.76517    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
```
