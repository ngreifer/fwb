# Reproducibility and Parallelization with \`fwb\`

Reproducibility ensures re-running the same analysis yields identical
results. Because a random process is involved in generating the
bootstrap weights with
[`fwb::fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md), steps
must be taken to ensure reproducibility is possible. As of version
0.7.0, all that is required is to call
[`set.seed()`](https://rdrr.io/r/base/Random.html) before
[`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md), replacing
`{N}` with your favorite integer:

``` r

set.seed({N})

f.out <- fwb(.)
```

The same procedure works for
[`vcovFWB()`](https://ngreifer.github.io/fwb/reference/vcovFWB.md). You
do not need to select a particular random number generator, and you do
not need to do anything differently when parallelizing. The rest of this
guide explains what is and is not guaranteed, how to reproduce results
computed by earlier versions, and what to know about the parallel
backends.

### What the results depend on

The bootstrap weights are a function of the state of the random number
generator when
[`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) is called,
and nothing about *how* the work is carried out enters into them. None
of the following changes your results:

- `cl`, which controls whether and how parallelization takes place,
  whether it is unspecified, an integer, a `cluster` object, or
  `"future"`

- the number of workers, or the *future* plan in use

- `verbose`, i.e., whether a progress bar is displayed

- the random number generator in effect
  ([`RNGkind()`](https://rdrr.io/r/base/Random.html))

- whether the function supplied to `statistic` draws random numbers of
  its own

An analysis is therefore reproducible on a machine with a different
number of cores, or by a collaborator who prefers a different parallel
backend, and switching from `cl = 2` to `cl = "future"` to inspect a
slow bootstrap costs nothing. What does determine the results is the
seed together with the arguments that define the bootstrap itself: `R`,
`wtype`, `strata`, `cluster`, and `simple` (see [The `simple`
argument](#simple)). The single exception is described in [Using
*futurize*](#futurize).

### Why this works

Before any bootstrap replication runs,
[`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) reserves a
separate random number stream for each of the `R` replications, and each
replication draws its weights from the stream reserved for it. These are
L’Ecuyer-CMRG streams, generated with
[`parallel::nextRNGStream()`](https://rdrr.io/r/parallel/RngStream.html).
The streams are derived from the calling session’s generator state,
which is what [`set.seed()`](https://rdrr.io/r/base/Random.html) sets,
and they are recorded in the returned object.

Two useful things follow. The weights for a given replication depend on
the seed and on the replication number, and not on which worker happened
to run it or in what order, so a parallel backend cannot change them.
And because the seeds that generated the weights are recorded, the
weights that were used can be recovered exactly rather than inferred,
which is what allows
[`fwb.array()`](https://ngreifer.github.io/fwb/reference/fwb.array.md)
to return them without the backend that produced them (see [BCa
confidence intervals](#bca)).

### Reproducing results from earlier versions

Results from *fwb* 0.6.0 and earlier cannot be reproduced with version
0.7.0 or later when `simple = TRUE`, which is the default for every
`wtype` except `"multinom"`. To reproduce a specific earlier analysis,
install the version it was run under[^1]. Results with `simple = FALSE`
are unchanged.

Previously, weights were drawn inside each bootstrap replication, from
whatever random number stream the worker running it happened to have.
This made the results depend on how the work was distributed, and
reproducing them required arranging the workers’ streams by hand. The
guidance was correspondingly involved, as it depended on `cl`, on
`simple`, on whether `statistic` was random, and on `verbose`, and with
a `cluster` object it required
[`parallel::clusterSetRNGStream()`](https://rdrr.io/r/parallel/RngStream.html)
instead of [`set.seed()`](https://rdrr.io/r/base/Random.html).

Two pieces of that guidance are worth calling out, because following
them today would be counterproductive.
[`parallel::clusterSetRNGStream()`](https://rdrr.io/r/parallel/RngStream.html)
no longer has any effect on
[`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md); the workers’
own streams are never used to draw weights, so code that relies on it in
place of [`set.seed()`](https://rdrr.io/r/base/Random.html) is not
reproducible. And `set.seed(., "L'Ecuyer-CMRG")` is no longer needed, as
[`set.seed()`](https://rdrr.io/r/base/Random.html) with any generator is
enough. Note that different generators select different streams, so the
same integer gives different (individually reproducible) results with
and without the `kind` argument.

### The `simple` argument

`simple` controls when the weights are drawn. With `simple = FALSE`, the
entire `R` \\\times\\ \\n\\ weight matrix is drawn in one step before
`statistic` is applied to anything. With `simple = TRUE`, each
replication draws its own weights as it runs, which avoids holding the
whole matrix in memory.

Both settings are fully reproducible, and both draw from the same weight
distribution, so either supports valid inference. They consume their
random number streams differently, though, so `simple = TRUE` and
`simple = FALSE` do not yield the same replicates as each other even
with the same seed. Previously they did coincide for the continuous
weight types, which is worth knowing if you have code that switched
between them freely. Default is `TRUE` except when `wtype = "multinom"`,
where `FALSE` is used so that results match `boot::boot(., stype = "f")`
exactly.

### BCa confidence intervals

Bias-corrected and accelerated (BCa) intervals require the empirical
influence of each unit on the bootstrap estimates, which means
recovering the weights that were used. Because those weights are
recorded (i.e., through the seeds that generated them),
`fwb.ci(., type = "bca")` and `summary(., ci.type = "bca")` work without
any preparation on your part, whatever `cl` was set to and whether or
not `statistic` is random. Previously, BCa intervals were refused
outright when `simple = TRUE` and `statistic` had a random component.

BCa intervals still require more bootstrap replications than there are
units in the data, and still cannot be used with `cluster` (i.e., with
the clustered bootstrap).

### Choosing a parallel backend

Because the choice does not affect your results, it can be made on
practical grounds. `cl` accepts an integer (which uses forking and is
not available on Windows), a `cluster` object from
[`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html),
or the string `"future"` to use whatever backend
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
has set up. The *futurize* package offers a further route, described
below.

One thing does differ between backends, and it does not concern
reproducibility. A `cluster` object starts fresh R processes with no
packages attached, so a `statistic` that calls functions from another
package needs that package attached there, e.g., with
`parallel::clusterEvalQ(cl, library(mypkg))`. *fwb* itself is attached
for you, so the `w_*()` functions (e.g.,
[`w_mean()`](https://ngreifer.github.io/fwb/reference/w_mean.md),
[`w_std()`](https://ngreifer.github.io/fwb/reference/w_mean.md)) can be
used with no extra effort. Forking and *future* backends handle this on
their own. Progress reporting also differs; see [Progress
bars](#progress).

### Using *futurize*

The *futurize* package can parallelize an
[`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) call by
piping it. See below for a code example:

``` r

library(futurize)
plan(multisession)

set.seed({N})

f.out <- fwb(.) |> futurize()
```

This is reproducible in the ordinary way, and a plain
[`set.seed()`](https://rdrr.io/r/base/Random.html) is all it needs.
There is one wrinkle, and it is the single exception to the rule that
only the seed matters: adding or removing `|> futurize()` changes the
replicates. *futurize* switches the generator to L’Ecuyer-CMRG for the
duration of the call, and doing so re-initializes the generator state
even when L’Ecuyer-CMRG is already in use. The streams *fwb* then
reserves are a deterministic function of your seed, but they are not the
ones a call without the pipe would have used, so reproducing a result
means reproducing whether
[`futurize()`](https://futurize.futureverse.org/reference/futurize.html)
was used.

See
[`vignette("futurize-81-fwb", package = "futurize")`](https://futurize.futureverse.org/articles/futurize-81-fwb.html)
for more details on using *futurize* with *fwb*.

### Progress bars

[`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md) displays a
progress bar when `verbose = TRUE`, using the *pbapply* package. This
works with every backend and, unlike in earlier versions, has no effect
on the results. Default is `TRUE` when running sequentially and `FALSE`
when parallelizing.

Alternatively, the *progressify* package routes progress through
*progressr*, which lets you choose the reporter. See below for a code
example:

``` r

library(progressify)
handlers(global = TRUE)

set.seed({N})

f.out <- fwb(.) |> progressify()
```

When combining parallelization with *progressify*, a *future* backend
must be used (i.e., `cl = "future"`, or the
[`futurize()`](https://futurize.futureverse.org/reference/futurize.html)
pipe described above), because *progressr* has no way to receive updates
from a forked worker. With `cl` set to an integer or a `cluster` object,
the bootstrap runs correctly but no progress is reported; use
`verbose = TRUE` in that case to display a progress bar.

See
[`vignette("progressify-81-fwb", package = "progressify")`](https://progressify.futureverse.org/articles/progressify-81-fwb.html)
for more details on using *progressify* with *fwb*.

### A note on `statistic`

`statistic` may draw random numbers, which no longer costs anything in
reproducibility, because each replication’s weights are drawn from a
stream reserved in advance and cannot be shifted by what `statistic` did
in an earlier replication. Its own draws are reproducible for the same
reason.

Two things are still worth avoiding. Do not call
[`set.seed()`](https://rdrr.io/r/base/Random.html) inside `statistic`,
which would give every replication the same draws. And note that a
`statistic` that parallelizes internally, on top of *fwb*’s own
parallelization, may not be reproducible; that is a property of the
nested code and not of
[`fwb()`](https://ngreifer.github.io/fwb/reference/fwb.md).

[^1]: To install a prior version of *fwb*, you can use, e.g.,
    `pak::pak("fwb@0.6.0")`.
