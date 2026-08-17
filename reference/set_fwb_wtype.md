# Set weights type

Set the default for the type of weights used in the weighted bootstrap
computed by [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md)
and [`vcovFWB()`](https://ngreifer.github.io/fwb/reference/vcovFWB.md).

## Usage

``` r
set_fwb_wtype(wtype = getOption("fwb_wtype", "exp"))

get_fwb_wtype(fwb)
```

## Arguments

- wtype:

  string; the type of weights to use. Allowable options include `"exp"`
  (the default), `"poisson"`, `"multinom"`, `"mammen"`, `"beta"`, and
  `"power"`. Abbreviations allowed. See
  [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) for what
  these mean.

- fwb:

  optional; an `<fwb>` object, the output of a call to
  [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md). If left
  empty, will extract the weights type from
  [`options()`](https://rdrr.io/r/base/options.html).

## Value

`set_fwb_wtype()` returns a call to
[`options()`](https://rdrr.io/r/base/options.html) with the options set
to those prior to `set_fwb_wtype()` being called. This makes it so that
calling `options(op)`, where `op` is the output of `set_fwb_wtype()`,
resets the `fwb_wtype` to its original value. `get_fwb_wtype()` returns
a string containing the `fwb_wtype` value set globally (if no argument
is supplied) or used in the supplied `<fwb>` object.

## Details

`set_fwb_wtype(x)` is equivalent to calling `options(fwb_wtype = x)`.
`get_fwb_wtype()` is equivalent to calling `getOption("fwb_wtype")` when
no argument is supplied and to extracting the `wtype` component of an
`<fwb>` object when supplied.

## See also

[`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) for a
definition of each types of weights;
[`vcovFWB()`](https://ngreifer.github.io/fwb/reference/vcovFWB.md);
[`options()`](https://rdrr.io/r/base/options.html);
[`boot::boot()`](https://rdrr.io/pkg/boot/man/boot.html) for the
traditional bootstrap.

## Examples

``` r
# Performing a Weibull analysis of the Bearing Cage
# failure data as done in Xu et al. (2020)
set.seed(123)
data("bearingcage")

#Set fwb type to "beta"
op <- set_fwb_wtype("beta")

weibull_est <- function(data, w) {
  fit <- survival::survreg(survival::Surv(hours, failure) ~ 1,
                           data = data, weights = w,
                           dist = "weibull")

  c(eta = unname(exp(coef(fit))), beta = 1/fit$scale)
}

boot_est <- fwb(bearingcage, statistic = weibull_est,
                R = 199, verbose = FALSE)
boot_est
#> FRACTIONAL WEIGHTED BOOTSTRAP
#> 
#> Call:
#> fwb(data = bearingcage, statistic = weibull_est, R = 199, verbose = FALSE)
#> 
#> Bootstrap Statistics :
#>          original         bias   std. error
#> eta  11792.178173 14292.864684 53929.165455
#> beta     2.035319     0.314054     1.144012

#Get the fwb type used in the bootstrap
get_fwb_wtype(boot_est)
#> [1] "beta"
get_fwb_wtype()
#> [1] "beta"

#Restore original options
options(op)

get_fwb_wtype()
#> [1] "exp"
```
