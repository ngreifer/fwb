#Integration with *progressify*, which rewrites `fwb(...)` so that progress is
#reported through *progressr* instead of *pbapply*.
#
#The transpiler reaches deep into `fwb()`'s calling convention: it wraps
#`statistic` in a function that signals one update per call, forces
#`verbose = FALSE`, passes the original statistic and the progressor through
#`fwb()`'s `...`, and -- crucially -- decides how many of the leading `statistic`
#calls to *skip* from `simple` and `wtype`, because `fwb()` calls `statistic` on
#unit weights once before the bootstrap (twice when `simple = TRUE`, to detect
#randomness).
#
#That skip count is the fragile part. If `fwb()` ever changes how many calibration
#calls it makes, or how `simple` defaults from `wtype`, the progress bar silently
#miscounts. The update-count assertions below are what would catch it.

skip_if_no_progressify <- function() {
  skip_on_cran()
  skip_if_not_installed("progressify")
  skip_if_not_installed("progressr")

  invisible(NULL)
}

test_that("`progressify()` supports `fwb()` and leaves the estimates unchanged", {
  skip_if_no_progressify()

  d <- fwb_test_df()

  expect_true("fwb" %in% progressify::progressify_supported_functions("fwb"))

  set.seed(2, "L'Ecuyer-CMRG")
  ref <- fwb(d, lm_stat, R = 20, verbose = FALSE)

  set.seed(2, "L'Ecuyer-CMRG")
  expect_no_error({
    p <- fwb(d, lm_stat, R = 20) |> progressify::progressify()
  })

  #The whole point of the wrapper is that it is numerically transparent: it must
  #not consume a draw or reorder anything.
  expect_s3_class(p, "fwb")
  expect_same_t(p, ref)
  expect_equal(p[["t0"]], ref[["t0"]], tolerance = fwb_eps())
  expect_identical(p[["R"]], ref[["R"]])
  expect_identical(dim(p[["t"]]), dim(ref[["t"]]))
  expect_identical(colnames(p[["t"]]), colnames(ref[["t"]]))
})

test_that("`progressify()` signals exactly `R` updates", {
  skip_if_no_progressify()

  d <- fwb_test_df()

  #`simple = FALSE`: one unit-weight calibration call, which the wrapper swallows.
  set.seed(2, "L'Ecuyer-CMRG")
  n <- count_progress_updates({
    invisible(fwb(d, lm_stat, R = 20, simple = FALSE) |> progressify::progressify())
  })

  expect_identical(n, 20L)

  #`wtype = "multinom"` flips the `simple` default to `FALSE`. The transpiler
  #re-derives the skip count from `wtype` when `simple` is absent, so this exercises
  #a third path through the same logic.
  set.seed(2, "L'Ecuyer-CMRG")
  n <- count_progress_updates({
    invisible(fwb(d, lm_stat, R = 20, wtype = "multinom") |>
                progressify::progressify())
  })

  expect_identical(n, 20L)

  #`R` supplied positionally rather than by name: the transpiler locates it by
  #matching against `fwb()`'s formals, which only works while `R` stays the third
  #argument.
  set.seed(2, "L'Ecuyer-CMRG")
  n <- count_progress_updates({
    invisible(fwb(d, lm_stat, 15, simple = FALSE) |> progressify::progressify())
  })

  expect_identical(n, 15L)
})

test_that("`progressify()` counts every replicate when `simple = TRUE`", {
  #*progressify*'s transpiler decides how many leading `statistic` calls to swallow from
  #`simple` and `wtype`:
  #
  #    .progressr_skip_count <- if (.progressr_simple) 2L else 1L
  #
  #The `2L` matched `fwb()` before version 0.7.0, which called `statistic` on unit
  #weights twice when `simple = TRUE` -- once for `t0` and once more to detect whether
  #`statistic` was random. Recording a stream per replicate made that detection
  #unnecessary, so there is now exactly one unit-weight call whatever `simple` is, and
  #the wrapper swallows one bootstrap replicate's update along with it.
  #
  #The effect is one missing tick out of `R`, so it is cosmetic -- but it is ours to
  #report: the fix is `1L` unconditionally, in *progressify*'s `fwb` transpiler.
  skip("upstream: progressify's skip count assumes two unit-weight `statistic` calls; report at github.com/futureverse/progressify")

  skip_if_no_progressify()

  d <- fwb_test_df()

  set.seed(2, "L'Ecuyer-CMRG")
  n <- count_progress_updates({
    invisible(fwb(d, lm_stat, R = 20) |> progressify::progressify())
  })

  expect_identical(n, 20L)
})

test_that("`progressify()` turns off the *pbapply* progress bar", {
  skip_if_no_progressify()

  d <- fwb_test_df()

  #Two reporters at once would interleave into nonsense, so the transpiler sets
  #`verbose = FALSE`. `fwb()`'s own default here would be `TRUE` (sequential).
  set.seed(2, "L'Ecuyer-CMRG")
  out <- utils::capture.output({
    invisible(count_progress_updates({
      invisible(fwb(d, lm_stat, R = 12) |> progressify::progressify())
    }))
  }, type = "output")

  expect_false(any(nzchar(out)))

  #An explicit `verbose = TRUE` is respected rather than overwritten -- the
  #transpiler only fills in the default.
  set.seed(2, "L'Ecuyer-CMRG")
  out <- utils::capture.output({
    invisible(fwb(d, lm_stat, R = 12, verbose = TRUE) |>
                progressify::progressify())
  }, type = "output")

  expect_true(any(nzchar(out)))
})

test_that("`progressify()` reports progress over a `future` backend", {
  skip_if_no_progressify()
  skip_if_not_installed("future.apply")

  d <- fwb_test_df()

  local_plan("multisession", workers = 2L)

  #*progressr* relays updates from `future` workers back to the parent, so the
  #count must survive the round trip.
  set.seed(2, "L'Ecuyer-CMRG")
  n <- count_progress_updates({
    invisible(fwb(d, lm_stat, R = 20, cl = "future", simple = FALSE) |>
                progressify::progressify())
  })

  expect_identical(n, 20L)

  set.seed(2, "L'Ecuyer-CMRG")
  ref <- fwb(d, lm_stat, R = 20, verbose = FALSE, cl = "future")

  set.seed(2, "L'Ecuyer-CMRG")
  p <- fwb(d, lm_stat, R = 20, cl = "future") |> progressify::progressify()

  expect_same_t(p, ref)
})

test_that("`progressify()` composes with `futurize()` in either order", {
  skip_if_no_progressify()
  skip_if_not_installed("futurize")
  skip_if_not_installed("future.apply")

  d <- fwb_test_df()

  local_plan("multisession", workers = 2L)

  set.seed(2, "L'Ecuyer-CMRG")
  n1 <- count_progress_updates({
    invisible(fwb(d, lm_stat, R = 20, simple = FALSE) |>
                progressify::progressify() |>
                futurize::futurize())
  })

  set.seed(2, "L'Ecuyer-CMRG")
  n2 <- count_progress_updates({
    invisible(fwb(d, lm_stat, R = 20, simple = FALSE) |>
                futurize::futurize() |>
                progressify::progressify())
  })

  expect_identical(n1, 20L)
  expect_identical(n2, 20L)

  #Adding progress reporting must not change the numbers, in either pipe order.
  set.seed(2, "L'Ecuyer-CMRG")
  ref <- fwb(d, lm_stat, R = 20, verbose = FALSE) |> futurize::futurize()

  set.seed(2, "L'Ecuyer-CMRG")
  a <- fwb(d, lm_stat, R = 20) |>
    progressify::progressify() |>
    futurize::futurize()

  set.seed(2, "L'Ecuyer-CMRG")
  b <- fwb(d, lm_stat, R = 20) |>
    futurize::futurize() |>
    progressify::progressify()

  expect_same_t(a, ref)
  expect_same_t(b, ref)
})

test_that("`progressify()` is reproducible", {
  skip_if_no_progressify()

  d <- fwb_test_df()

  set.seed(17, "L'Ecuyer-CMRG")
  a <- fwb(d, lm_stat, R = 20) |> progressify::progressify()

  set.seed(17, "L'Ecuyer-CMRG")
  b <- fwb(d, lm_stat, R = 20) |> progressify::progressify()

  expect_same_t(a, b)

  set.seed(18, "L'Ecuyer-CMRG")
  other <- fwb(d, lm_stat, R = 20) |> progressify::progressify()

  expect_different_t(other, a)
})

test_that("`progressify()` passes extra `statistic` arguments through", {
  skip_if_no_progressify()

  d <- fwb_test_df()

  #The transpiler smuggles `...FUN` and `.progressr_progressor` through `fwb()`'s
  #`...`, which is also how the user's own extra arguments get to `statistic`. If
  #`fwb()` ever started validating `...` against `statistic`'s formals, this is
  #where it would break.
  scaled_stat <- function(data, w, mult = 1) {
    mult * coef(lm(y ~ x1 + x2, data = data, weights = w))
  }

  set.seed(19, "L'Ecuyer-CMRG")
  p <- fwb(d, scaled_stat, R = 12, mult = 3) |> progressify::progressify()

  set.seed(19, "L'Ecuyer-CMRG")
  ref <- fwb(d, scaled_stat, R = 12, verbose = FALSE, mult = 3)

  expect_same_t(p, ref)
  expect_equal(p[["t0"]], 3 * lm_stat(d, rep(1, nrow(d))),
               tolerance = fwb_eps(), ignore_attr = TRUE)
})

test_that("`progressify()` output supports the usual downstream methods", {
  skip_if_no_progressify()

  d <- fwb_test_df()

  #`R` is comfortably above `nrow(d)` so that the BCa interval `summary()` computes
  #does not land on the extreme order statistics and warn. That warning is real and is
  #tested in `test-fwb-ci.R`; here it would just be noise from an unrelated method.
  set.seed(20, "L'Ecuyer-CMRG")
  p <- fwb(d, lm_stat, R = 200) |> progressify::progressify()

  expect_output(print(p), "FRACTIONAL WEIGHTED BOOTSTRAP")

  expect_no_error({
    s <- summary(p)
  })

  expect_s3_class(s, "summary.fwb")

  #`fwb.array()` re-generates the weights from `$seed` and `$wtype` only, so it is
  #indifferent to the wrapped `statistic` -- but `$statistic` in the returned object
  #is now *progressify*'s wrapper rather than the user's function, which is worth
  #recording so that nothing downstream comes to rely on it being the original.
  expect_no_error({
    a <- fwb.array(p)
  })

  expect_identical(dim(a), c(200L, nrow(d)))

  expect_false(isTRUE(all.equal(p[["statistic"]], lm_stat)))
})

test_that("progress reporting for non-`future` backends is a known limitation", {
  skip_if_no_progressify()
  skip_if_no_forking()

  d <- fwb_test_df()

  #`cl = <integer>` forks, and *progressr* has no channel back from a forked
  #child, so no updates arrive. Recording it here means the day it starts working
  #is a visible change rather than a surprise -- and until then, `?fwb` should
  #point users at a `future` backend when they want a *progressr* bar.
  set.seed(2, "L'Ecuyer-CMRG")
  n <- count_progress_updates({
    invisible(fwb(d, lm_stat, R = 12, cl = 2) |> progressify::progressify())
  })

  expect_identical(n, 0L)

  #The estimates are still correct; only the reporting is missing.
  set.seed(2, "L'Ecuyer-CMRG")
  ref <- fwb(d, lm_stat, R = 12, verbose = FALSE, cl = 2)

  set.seed(2, "L'Ecuyer-CMRG")
  p <- fwb(d, lm_stat, R = 12, cl = 2) |> progressify::progressify()

  expect_same_t(p, ref)
})
