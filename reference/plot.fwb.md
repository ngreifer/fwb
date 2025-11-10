# Plots of the Output of a Fractional Weighted Bootstrap

`plot.fwb()` takes an `fwb` object and produces plots for the bootstrap
replicates of the statistic of interest.

## Usage

``` r
# S3 method for class 'fwb'
plot(
  x,
  index = 1L,
  qdist = "norm",
  nclass = NULL,
  df,
  type = c("hist", "qq"),
  ...
)
```

## Arguments

- x:

  an `fwb` object; the output of a call to
  [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md).

- index:

  the index of the position of the quantity of interest in `x$t0` if
  more than one was specified in
  [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md). Only one
  value is allowed at a time. By default the first statistic is used.

- qdist:

  `character`; when a Q-Q plot is requested (as it is by default; see
  `type` argument below), the distribution against which the Q-Q plot
  should be drawn. Allowable options include `"norm"` (normal
  distribution - the default) and `"chisq"` (chi-squared distribution).

- nclass:

  when a histogram is requested (as it is by default; see `type`
  argument below), the number of classes to be used. The default is the
  integer between 10 and 100 closest to `ceiling(length(R)/25)` where
  `R` is the number of bootstrap replicates.

- df:

  if `qdist` is `"chisq"`, the degrees of freedom for the chi-squared
  distribution to be used. If not supplied, the degrees of freedom will
  be estimated using maximum likelihood.

- type:

  the type of plot to display. Allowable options include `"hist"` for a
  histogram of the bootstrap estimates and `"qq"` for a Q-Q plot of the
  estimates against the distribution supplied to `qdist`.

- ...:

  ignored.

## Value

`x` is returned invisibly.

## Details

This function can produces two side-by-side plots: a histogram of the
bootstrap replicates and a Q-Q plot of the bootstrap replicates against
theoretical quantiles of a supplied distribution (normal or
chi-squared). For the histogram, a vertical dotted line indicates the
position of the estimate computed in the original sample. For the Q-Q
plot, the expected line is plotted.

## See also

[`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md),
[`summary.fwb()`](https://ngreifer.github.io/fwb/reference/summary.fwb.md),
[`boot::plot.boot()`](https://rdrr.io/pkg/boot/man/plot.boot.html),
[`hist()`](https://rdrr.io/r/graphics/hist.html),
[`qqplot()`](https://rdrr.io/r/stats/qqnorm.html)

## Examples

``` r
# See examples at help("fwb")
```
