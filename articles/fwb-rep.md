# Reproducibility and Parallelization with \`fwb\`

Reproducibility ensures re-running the same analysis yields identical
results. Because a random process is involved in generating the
bootstrap weights with
[`fwb::fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md), steps
must be taken to ensure reproducibility is possible.

There are a few arguments to
[`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) that are
relevant for reproducibility. These are `statistic`, `simple`, and `cl`.

- `statistic` is the function that is applied to each bootstrap dataset
  and returns the quantities of interest to be estimated. It can either
  have a random component or not; each case requires special attention.
  It is always safer to avoid having a random component in `statistic`.
  Most common regression functions to do not involve a random component,
  but some advanced models, like machine learning models, may. Do not
  ever include a call to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) or supply a seed to
  `statistic`. If you are using parallelization for
  [`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md), do not use
  parallelization within `statistic`.

- `simple` controls whether the bootstrap weights are generated all at
  once (`simple = FALSE`) or generated separately within each bootstrap
  iteration (`simple = TRUE`). When `simple = FALSE`, the weights are
  generated before any parallelization takes place or `statistic` is
  called, which makes ensuring reproducibility more straightforward.
  When `simple = TRUE`, the weights are generated before the call to
  `statistic` in each bootstrap iteration, which can make it a bit more
  challenging to ensure reproducibility when using parallelization and
  adds even more challenges when `statistic` also has a random
  component.

- `cl` controls whether and how parallelization takes place. It is
  passed directly to
  [`pbapply::pblapply()`](https://peter.solymos.org/pbapply/reference/pbapply.html),
  which calls either
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html),
  [`parallel::parLapply()`](https://rdrr.io/r/parallel/clusterApply.html),
  or
  [`future.apply::future_lapply()`](https://future.apply.futureverse.org/reference/future_lapply.html),
  depending on how it is specified. The usual arguments include an
  integer referring to the number of cores, which only works on Mac and
  triggers
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html); a
  `cluster` object (usually the result of a call to
  [`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html)
  or related functions), which triggers
  [`parallel::parLapply()`](https://rdrr.io/r/parallel/clusterApply.html);
  or `"future"`, which uses a `future` backend (usually initialized
  using
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)).
  Each of these involves different requirements for ensuring
  reproducibility.

This guide will proceed for combinations of these scenarios.

### Case 1: No parallelization (`cl = NULL`)

When no parallelization is used (i.e., `cl` is unspecified, `NULL`, or
`1`), all you need to do is call
[`set.seed()`](https://rdrr.io/r/base/Random.html) before
[`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) to ensure
reproducibility. It doesn’t matter what `simple` or `statistic` do. This
is probably the most common case. Just run the following to ensure
reproducibility, replacing `{N}` with your favorite integer.

``` r
set.seed({N})

f.out <- fwb(.)
```

### Case 2: `simple = FALSE`, non-random `statistic`

If `simple = FALSE` and `statistic` does not have a random component,
see Case 1, regardless of whether or how parallelization is used. In
this case, no random process occurs within each cluster, so no special
steps need to be taken beyond setting a seed. Note that `simple = TRUE`
by default unless `wtype = "multinom"`, so this must be set manually.
See below for a code example:

``` r
set.seed({N})

f.out <- fwb(., simple = FALSE)
```

### Case 3: `cl` is an integer

When `cl` is an integer and the criteria for Case 2 are not met (i.e.,
`simple = TRUE` or `statistic` has a random component), one additional
step is required for ensuring reproducibility. Again, all you need to do
is use [`set.seed()`](https://rdrr.io/r/base/Random.html), but you must
call it with `kind = "L'Ecuyer-CMRG"`, which is the only method
appropriate for use across multiple clusters. See below for a code
example:

``` r
set.seed({N}, "L'Ecuyer-CMRG")

f.out <- fwb(., cl = 3)
```

### Case 4: `cl` is `"future"`

When using a `future` backend and the criteria for Case are not met, you
can use the same solution as for Case 3.
[`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) performs an
additional step to make sure the seed is correctly sent to
[`future.apply::future_lapply()`](https://future.apply.futureverse.org/reference/future_lapply.html).
(Internally, this works by setting `future.seed = TRUE`, which you
should not do yourself.) See below for a code example:

``` r
library(future)

plan(multisession, workers = 3)
set.seed({N}, "L'Ecuyer-CMRG")

f.out <- fwb(., cl = "future")
```

### Case 5: `cl` is a `cluster` object

When `cl` is a `cluster` object (i.e., the output of a call to
[`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html),
[`parallel::makePSOCKcluster()`](https://rdrr.io/r/parallel/makeCluster.html),
[`parallel::makeForkCluster()`](https://rdrr.io/r/parallel/makeCluster.html)
or similar functions in *parallelly*), an additional step must be taken
to ensure reproducibility. Unfortunately, you can’t use
[`set.seed()`](https://rdrr.io/r/base/Random.html); you have to use
[`parallel::clusterSetRNGStream()`](https://rdrr.io/r/parallel/RngStream.html),
to which you supply the `cluster` object and your desired seed. See
below for a code example:

``` r
library(parallel)

cl <- makeCluster(3)
clusterSetRNGStream(cl, {N})

f.out <- fwb(., cl = cl)
```

## Computing BCa confidence intervals

Although the main purpose of considering reproducibility is to ensure
that multiple runs of the same code produce identical results, there is
another situation in which it can be important to be able to reproduce
the weights, and that is when computing bias-corrected accelerated (BCa)
confidence intervals using `fwb.ci(., type = "bca")` or
`summary(., ci.type = "bca")`. BCa confidence intervals have the best
statistical properties among the available bootstrap confidence
intervals, but they require computing the influence each unit has on the
bootstrap estimates, which requires re-generating the weights as they
were generated by
[`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md).

There are some cases where you don’t have to do any special work to
ensure BCa intervals are correctly computed. These include:

- `simple = FALSE`, regardless of parallelization or randomness in
  `statistic`
- `simple = TRUE`, there is no randomness in `statistic`, and no
  parallelization is used
- `simple = TRUE`, there is no randomness in `statistic`, and `cl` is an
  integer or `"future"`

In these cases,
[`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) saves the
state of the random seed that was used to originally generate the
weights, and
[`fwb.ci()`](https://ngreifer.github.io/fwb/reference/fwb.ci.md) recalls
that seed to re-generate the weights and then computes the required
statistics for the BCa interval without requiring any extra involvement
by the user. Remember, `simple = TRUE` by default unless
`wtype = "multinom"`.

Otherwise, when the following condition is met, an additional step is
required:

- `simple = TRUE`, there is no randomness in `statistic`, and `cl` is a
  `cluster` object

In this case, you need to call `parallel::clusterSetRNGStream(cl, {N})`
with the same seed as as was used prior to
[`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) immediately
before calling
[`fwb.ci()`](https://ngreifer.github.io/fwb/reference/fwb.ci.md) or
[`summary()`](https://rdrr.io/r/base/summary.html).

When `simple = TRUE` and there is any randomness in `statistic`, it is
not possible to re-generate the weights that were used in the bootstrap,
so BCa confidence intervals cannot be computed.
[`fwb.ci()`](https://ngreifer.github.io/fwb/reference/fwb.ci.md) (and
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`confint()`](https://rdrr.io/r/stats/confint.html), which call
[`fwb.ci()`](https://ngreifer.github.io/fwb/reference/fwb.ci.md))
automatically checks for this case and throws an error if BCa confidence
intervals are requested when these conditions are met.
