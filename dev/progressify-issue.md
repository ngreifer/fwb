## `fwb` transpiler: progress bar is one tick short when `simple = TRUE`

The `fwb` transpiler's skip count assumes `fwb()` calls `statistic` on unit weights
**twice** when `simple = TRUE`. In the current *fwb* development version it calls it once,
whatever `simple` is, so the wrapper swallows the first real bootstrap replicate's update
and the bar finishes one tick short of `R`.

I maintain *fwb* and made the change that broke this, so apologies for the breakage. This
is not yet on CRAN — the released 0.6.0 still behaves the way the transpiler expects, so
nothing is broken for users today. I'm reporting it ahead of the release rather than after.
Happy to coordinate on timing, and to add a test on my side once you've picked a fix.

### Reproducible example

```r
library(fwb)
library(progressify)

## Count progress updates without rendering a bar.
count_ticks <- function(expr) {
  n <- 0L
  progressr::with_progress(
    withCallingHandlers(force(expr),
                        condition = function(cnd) {
                          if (inherits(cnd, "progression") &&
                              identical(cnd[["type"]], "update")) n <<- n + 1L
                        }),
    handlers = progressr::handler_void(), enable = TRUE)
  n
}

set.seed(1)
d <- data.frame(x = rnorm(30))
d$y <- d$x + rnorm(30)
stat <- function(data, w) coef(lm(y ~ x, data = data, weights = w))

count_ticks(fwb(d, stat, R = 20, simple = TRUE)  |> progressify())  # 19, expected 20
count_ticks(fwb(d, stat, R = 20, simple = FALSE) |> progressify())  # 20, correct
count_ticks(fwb(d, stat, R = 20)                 |> progressify())  # 19 (simple defaults TRUE)
count_ticks(fwb(d, stat, R = 20, wtype = "multinom") |> progressify())  # 20 (simple defaults FALSE)
```

| call | expected | actual |
|---|---|---|
| `simple = TRUE` | 20 | **19** |
| `simple = FALSE` | 20 | 20 |
| default (`simple` unspecified) | 20 | **19** |
| `wtype = "multinom"` | 20 | 20 |

So it affects the default for every `wtype` except `"multinom"` — i.e. most calls.

### Cause

`template_fwb_outer` sets:

```r
.progressr_skip_count <- if (.progressr_simple) 2L else 1L
```

and `template_fwb_FUN` decrements that on each call before it starts ticking. The count is
the number of times `fwb()` calls `statistic` before the `R` replicates begin.

In *fwb* 0.6.0 and earlier, `simple = TRUE` did call `statistic` twice on unit weights: the
second call existed only to detect whether `statistic` had a random component, which
determined whether the bootstrap weights could later be recovered. The development version
recovers the weights from recorded RNG streams instead, so that detection became
unnecessary and the second call was removed. There is now exactly one unit-weight call in
both cases:

```r
for (simple in c(TRUE, FALSE)) {
  calls <- 0L
  counting_stat <- function(data, w) {
    if (all(w == 1)) calls <<- calls + 1L
    coef(lm(y ~ x, data = data, weights = w))
  }
  invisible(fwb(d, counting_stat, R = 20, simple = simple, verbose = FALSE))
  cat(sprintf("simple = %-5s  %d unit-weight call(s)\n", simple, calls))
}
#> simple = TRUE   1 unit-weight call(s)
#> simple = FALSE  1 unit-weight call(s)
```

### Suggested fix

An unconditional

```r
.progressr_skip_count <- 1L
```

is correct for the development version and for every release after it.

It would break *fwb* 0.6.0, though: `simple = TRUE` there would leave the second
unit-weight call to be counted as a tick, giving `R + 1` updates against a progressor with
`R` steps. So if you want to keep 0.6.0 working, this needs a version test:

```r
.progressr_skip_count <- if (.progressr_simple &&
                             utils::packageVersion("fwb") < "<next release>") 2L else 1L
```

I've deliberately left the version as a placeholder: the release number isn't settled on my
end yet, and I don't want to hand you one that turns out to be wrong. **Tell me which of
the two you'd prefer and I'll follow up with the exact version once it's tagged** — or if
you'd rather just require the new release and drop the branch entirely, that's fine by me
and I'll note the dependency in *fwb*'s docs.

Either way, note that an unconditional skip count also removes the need to infer `simple`,
which the template currently does by reproducing `fwb()`'s default rule:

```r
.progressr_simple <- if (is.null(SIMPLE)) WTYPE != "multinom" else SIMPLE
```

With a constant skip count, `SIMPLE` and `WTYPE` are not needed at all. That seems like a
robustness win independent of this bug: the rule is an implementation detail of `fwb()`
that the transpiler has no way to see changes to, and it has already drifted once.

### Related, separate: `R` is not always the number of replicates

Worth flagging while you're in this code, since it also makes `.progressr_steps` wrong.

The development version adds an exhaustive multinomial bootstrap: because there are only
`n^n` distinct multinomial resamples of `n` units, when `R` reaches that many `fwb()` uses
each exactly once and **reduces `R`** to that number. The transpiler takes `steps` from the
`R` in the call, so the bar stalls partway:

```r
d2 <- data.frame(x = c(1, 2, 4), y = c(1, 2, 3))
stat2 <- function(data, w) c(m = w_mean(data$x, w))
count_ticks(fwb(d2, stat2, R = 999, wtype = "multinom") |> progressify())
#> 27   -- progressor was created with 999 steps
```

Only reachable for very small `n` (`n = 6` already has 46656 resamples, so with the default
`R = 999` nothing above `n = 4` is affected), and `fwb()` emits a message saying it has done
this. There may be nothing worth doing about it — the number of replicates isn't knowable
from the call alone, so a transpiler arguably cannot get this right in general. Recording it
rather than proposing a fix.

### Session info

```
fwb 0.6.0.9001 (development, github.com/ngreifer/fwb)
progressify 0.2.0
progressr 1.0.0
R 4.6.1, macOS
```
