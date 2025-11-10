# Fractional Weighted Bootstrap Covariance Matrix Estimation

`vcovFWB()` estimates the covariance matrix of model coefficient
estimates using the fractional weighted bootstrap. It serves as a
drop-in for `[stats::vcov()]` or
[`sandwich::vcovBS()`](https://sandwich.R-Forge.R-project.org/reference/vcovBS.html).
Clustered covariances are can be requested.

## Usage

``` r
vcovFWB(
  x,
  cluster = NULL,
  R = 1000,
  start = FALSE,
  wtype = getOption("fwb_wtype", "exp"),
  drop0 = FALSE,
  ...,
  fix = FALSE,
  use = "pairwise.complete.obs",
  .coef = stats::coef,
  verbose = FALSE,
  cl = NULL
)
```

## Arguments

- x:

  a fitted model object, such as the output of a call to
  [`lm()`](https://rdrr.io/r/stats/lm.html) or
  [`glm()`](https://rdrr.io/r/stats/glm.html). The model object must
  result from a function that can be updated using
  [`update()`](https://rdrr.io/r/stats/update.html) and has a `weights`
  argument to input non-integer case weights.

- cluster:

  a variable indicating the clustering of observations, a `list` (or
  `data.frame`) thereof, or a formula specifying which variables from
  the fitted model should be used (see examples). By default
  (`cluster = NULL`), either `attr(x, "cluster")` is used (if any) or
  otherwise every observation is assumed to be its own cluster.

- R:

  the number of bootstrap replications. Default is 1000 (more is better
  but slower).

- start:

  `logical`; should `.coef(x)` be passed as `start` to the
  `update(x, weights = ...)` call? In case the model `x` is computed by
  some numeric iteration, this may speed up the bootstrapping. Default
  is `FALSE`.

- wtype:

  string; the type of weights to use. Allowable options include `"exp"`
  (the default), `"pois"`, `"multinom"`, `"mammen"`, `"beta"`, and
  `"power"`. See
  [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) for
  details. See
  [`set_fwb_wtype()`](https://ngreifer.github.io/fwb/reference/set_fwb_wtype.md)
  to set a global default.

- drop0:

  `logical`; when `wtype` is `"multinom"` or `"poisson"`, whether to
  drop units that are given weights of 0 from the model call in each
  iteration. If `TRUE`, the model will be called with an additional
  `subset` argument, filtering out units with weights of 0 (note this
  will overwrite any argument to `subset` in the original call). If
  `NA`, weights of 0 will be set to `NA` instead. Ignored for other
  `wtype`s because they don't produce 0 weights. Default is `FALSE`.

- ...:

  ignored.

- fix:

  `logical`; if `TRUE`, the covariance matrix is fixed to be positive
  semi-definite in case it is not.

- use:

  `character`; specification passed to
  [`stats::cov()`](https://rdrr.io/r/stats/cor.html) for handling
  missing coefficients/parameters.

- .coef:

  a function used to extract the coefficients from each fitted model.
  Must return a numeric vector. By default,
  [`stats::coef`](https://rdrr.io/r/stats/coef.html) is used, but
  `marginaleffects::get_coef` can be a more reliable choice for some
  models that have a non-standard
  [`coef()`](https://rdrr.io/r/stats/coef.html) method, like that for
  [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html)
  models.

- verbose:

  `logical`; whether to display a progress bar.

- cl:

  a cluster object created by
  [`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html),
  an integer to indicate the number of child-processes (integer values
  are ignored on Windows) for parallel evaluations, or the string
  `"future"` to use a `future` backend. See the `cl` argument of
  [`pbapply::pblapply()`](https://peter.solymos.org/pbapply/reference/pbapply.html)
  for details. If `NULL`, no parallelization will take place. See
  [`vignette("fwb-rep")`](https://ngreifer.github.io/fwb/articles/fwb-rep.md)
  for details.

## Value

A matrix containing the covariance matrix estimate.

## Details

`vcovFWB()` functions like other
[`vcov()`](https://rdrr.io/r/stats/vcov.html)-like functions, such as
those in the sandwich package, in particular,
[`sandwich::vcovBS()`](https://sandwich.R-Forge.R-project.org/reference/vcovBS.html),
which implements the traditional bootstrap (and a few other bootstrap
varieties for linear models). Sets of weights are generated as described
in the documentation for
[`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md), and the
supplied model is re-fit using those weights. When the fitted model
already has weights, these are multiplied by the bootstrap weights.

For `lm` objects, the model is re-fit using
[`.lm.fit()`](https://rdrr.io/r/stats/lmfit.html) for speed, and,
similarly, `glm` objects are re-fit using
[`glm.fit()`](https://rdrr.io/r/stats/glm.html) (or whichever fitting
method was originally used). For other objects,
[`update()`](https://rdrr.io/r/stats/update.html) is used to populate
the weights and re-fit the model (this assumes the fitting function
accepts non-integer case weights through a `weights` argument). If a
model accepts weights in some other way,
[`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) should be
used instead; `vcovFWB()` is inherently limited in its ability to handle
all possible models. It is important that the original model was not fit
using frequency weights (i.e., weights that allow one row of data to
represent multiple full, identical, individual units) unless clustering
is used.

See
[`sandwich::vcovBS()`](https://sandwich.R-Forge.R-project.org/reference/vcovBS.html)
and
[`sandwich::vcovCL()`](https://sandwich.R-Forge.R-project.org/reference/vcovCL.html)
for more information on clustering covariance matrices, and see
[`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) for more
information on how clusters work with the fractional weighted bootstrap.
When clusters are specified, each cluster is given a bootstrap weight,
and all members of the cluster are given that weight; estimation then
proceeds as normal. By default, when `cluster` is unspecified, each unit
is considered its own cluster.

## See also

- [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) for
  performing the fractional weighted bootstrap on an arbitrary quantity

- [`fwb.ci()`](https://ngreifer.github.io/fwb/reference/fwb.ci.md) for
  computing nonparametric confidence intervals for `fwb` objects

- [`summary.fwb()`](https://ngreifer.github.io/fwb/reference/summary.fwb.md)
  for producing standard errors and confidence intervals for `fwb`
  objects

- [`sandwich::vcovBS()`](https://sandwich.R-Forge.R-project.org/reference/vcovBS.html)
  for computing covariance matrices using the traditional bootstrap (the
  fractional weighted bootstrap is also available but with limited
  options).

## Examples

``` r
set.seed(123, "L'Ecuyer-CMRG")
data("infert")
fit <- glm(case ~ spontaneous + induced, data = infert,
             family = "binomial")
lmtest::coeftest(fit, vcov. = vcovFWB, R = 200)
#> Warning: non-integer #successes in a binomial glm!
#> 
#> z test of coefficients:
#> 
#>             Estimate Std. Error z value  Pr(>|z|)    
#> (Intercept) -1.70786    0.26445 -6.4581 1.060e-10 ***
#> spontaneous  1.19721    0.21868  5.4746 4.384e-08 ***
#> induced      0.41813    0.19723  2.1201     0.034 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
# Example from help("vcovBS", package = "sandwich")
data("PetersenCL", package = "sandwich")
m <- lm(y ~ x, data = PetersenCL)

# Note: this is not to compare performance, just to
# demonstrate the syntax
cbind(
  "BS" = sqrt(diag(sandwich::vcovBS(m))),
  "FWB" = sqrt(diag(vcovFWB(m))),
  "BS-cluster" = sqrt(diag(sandwich::vcovBS(m, cluster = ~firm))),
  "FWB-cluster" = sqrt(diag(vcovFWB(m, cluster = ~firm)))
)
#>                     BS        FWB BS-cluster FWB-cluster
#> (Intercept) 0.02834177 0.02862498 0.07244951  0.06554226
#> x           0.02642269 0.02837431 0.04724351  0.04789381

# Using `wtype = "multinom"` exactly reproduces
# `sandwich::vcovBS()`
set.seed(11)
s <- sandwich::vcovBS(m, R = 200)
set.seed(11)
f <- vcovFWB(m, R = 200, wtype = "multinom")

all.equal(s, f)
#> [1] TRUE
# Using a custom argument to `.coef`
set.seed(123)
data("infert")

fit <- nnet::multinom(education ~ age, data = infert,
                      trace = FALSE)

# vcovFWB(fit, R = 200) ## error
coef(fit) # coef() returns a matrix
#>         (Intercept)         age
#> 6-11yrs    5.482932 -0.09335993
#> 12+ yrs    9.372401 -0.21892472

# Write a custom function to extract vector of
# coefficients (can also use marginaleffects::get_coef())
coef_multinom <- function(x) {
  p <- t(coef(x))

  setNames(as.vector(p),
           paste(colnames(p)[col(p)],
                 rownames(p)[row(p)],
                 sep = ":"))
}
coef_multinom(fit) # returns a vector
#> 6-11yrs:(Intercept)         6-11yrs:age 12+ yrs:(Intercept)         12+ yrs:age 
#>          5.48293157         -0.09335993          9.37240097         -0.21892472 

vcovFWB(fit, R = 200, .coef = coef_multinom)
#>                     6-11yrs:(Intercept)  6-11yrs:age 12+ yrs:(Intercept)
#> 6-11yrs:(Intercept)           7.1290604 -0.202821220           6.7927253
#> 6-11yrs:age                  -0.2028212  0.005837787          -0.1923858
#> 12+ yrs:(Intercept)           6.7927253 -0.192385799           7.1922051
#> 12+ yrs:age                  -0.1922646  0.005504376          -0.2046696
#>                      12+ yrs:age
#> 6-11yrs:(Intercept) -0.192264594
#> 6-11yrs:age          0.005504376
#> 12+ yrs:(Intercept) -0.204669617
#> 12+ yrs:age          0.005897351
```
