# Extract Confidence Intervals from a `bootci` Object

`get_ci()` extracts the confidence intervals from the output of a call
to [`boot::boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html) or
[`fwb.ci()`](https://ngreifer.github.io/fwb/reference/fwb.ci.md) in a
clean way. Normally the confidence intervals can be a bit challenging to
extract because of the unusual structure of the object.

## Usage

``` r
get_ci(x, type = "all")
```

## Arguments

- x:

  a `bootci` object; the output of a call to
  [`boot::boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html) or
  [`fwb.ci()`](https://ngreifer.github.io/fwb/reference/fwb.ci.md).

- type:

  the type of confidence intervals to extract. Only those available in
  `x` are allowed. Should be a given as a subset of the types passed to
  `type` in `boot.ci()` or
  [`fwb.ci()`](https://ngreifer.github.io/fwb/reference/fwb.ci.md). The
  default, `"all"`, extracts all confidence intervals in `x`.

## Value

A list with an entry for each confidence interval type; each entry is a
numeric vector of length 2 with names `"L"` and `"U"` for the lower and
upper interval bounds, respectively. The `"conf"` attribute contains the
confidence level.

## See also

[`fwb.ci()`](https://ngreifer.github.io/fwb/reference/fwb.ci.md),
[`confint.fwb()`](https://ngreifer.github.io/fwb/reference/summary.fwb.md),
[`boot::boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html)

## Examples

``` r
#See example at help("fwb.ci")
```
