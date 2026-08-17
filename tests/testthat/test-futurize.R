#Integration with *futurize*, which rewrites `fwb(...)` into
#`fwb(..., cl = "future")` when a call is piped to `futurize()`.
#
#The transpiler does three things that matter to *fwb*: it appends `cl = "future"`,
#it sets `future.seed = TRUE` through `pbapply`, and it temporarily switches the
#RNG to L'Ecuyer-CMRG for the duration of the call. That last step is the reason
#`futurize()` is reproducible even from a Mersenne-Twister stream -- and also the
#reason its replicates do not match a hand-written `cl = "future"` call.
#
#These tests exist because that contract is enforced entirely from *futurize*'s
#side: nothing in *fwb* would notice if `fwb()`'s formals or its use of `pbapply`
#drifted out from under it.

skip_if_no_futurize <- function() {
  skip_on_cran()
  skip_if_not_installed("futurize")
  skip_if_not_installed("future.apply")

  invisible(NULL)
}

test_that("`futurize()` transpiles `fwb()` to `cl = \"future\"`", {
  skip_if_no_futurize()

  d <- fwb_test_df()

  #`fwb()` and `vcovFWB()` are the two functions *futurize* claims to support; if
  #either is dropped, the vignette cross-reference in `?fwb` goes stale.
  expect_true(all(c("fwb", "vcovFWB") %in%
                    futurize::futurize_supported_functions("fwb")))

  local_plan("multisession", workers = 2L)

  set.seed(4, "L'Ecuyer-CMRG")
  expect_no_error({
    f <- fwb(d, lm_stat, R = 20, verbose = FALSE) |> futurize::futurize()
  })

  expect_s3_class(f, "fwb")
  expect_identical(f[["call"]][["cl"]], "future")

  #The backend is not recorded in the object: weights come from streams recorded before
  #dispatch, so nothing downstream needs to know how the run was parallelized.
  expect_null(attr(f, "cl", TRUE))
})

test_that("`futurize()` gives reproducible results", {
  skip_if_no_futurize()

  d <- fwb_test_df()

  local_plan("multisession", workers = 2L)

  set.seed(4, "L'Ecuyer-CMRG")
  a <- fwb(d, lm_stat, R = 20, verbose = FALSE) |> futurize::futurize()

  set.seed(4, "L'Ecuyer-CMRG")
  b <- fwb(d, lm_stat, R = 20, verbose = FALSE) |> futurize::futurize()

  expect_same_t(a, b)

  set.seed(5, "L'Ecuyer-CMRG")
  other <- fwb(d, lm_stat, R = 20, verbose = FALSE) |> futurize::futurize()

  expect_different_t(other, a)

  #*futurize* switches to L'Ecuyer-CMRG itself, so a plain `set.seed()` is enough
  #for reproducibility -- unlike a hand-written `cl = "future"`, which needs the
  #`kind` argument. Worth pinning: it is the main ergonomic benefit of the pipe.
  old_kind <- RNGkind("Mersenne-Twister")
  defer_call(quote(suppressWarnings(RNGkind(.))), old_kind[[1L]])

  set.seed(6)
  m1 <- fwb(d, lm_stat, R = 20, verbose = FALSE) |> futurize::futurize()

  RNGkind("Mersenne-Twister")
  set.seed(6)
  m2 <- fwb(d, lm_stat, R = 20, verbose = FALSE) |> futurize::futurize()

  expect_same_t(m2, m1)

  RNGkind("Mersenne-Twister")
  set.seed(7)
  m3 <- fwb(d, lm_stat, R = 20, verbose = FALSE) |> futurize::futurize()

  expect_different_t(m3, m1)
})

test_that("`futurize()` results do not depend on the plan or worker count", {
  skip_if_no_futurize()

  d <- fwb_test_df()

  local_plan("multisession", workers = 2L)

  set.seed(4, "L'Ecuyer-CMRG")
  w2 <- fwb(d, lm_stat, R = 20, verbose = FALSE) |> futurize::futurize()

  local_plan("multisession", workers = max_test_workers())

  set.seed(4, "L'Ecuyer-CMRG")
  w3 <- fwb(d, lm_stat, R = 20, verbose = FALSE) |> futurize::futurize()

  expect_same_t(w3, w2)

  local_plan("sequential")

  set.seed(4, "L'Ecuyer-CMRG")
  s <- fwb(d, lm_stat, R = 20, verbose = FALSE) |> futurize::futurize()

  expect_same_t(s, w2)
})

test_that("`futurize()` restores the RNG kind afterwards", {
  skip_if_no_futurize()

  d <- fwb_test_df()

  local_plan("sequential")

  old_kind <- RNGkind("Mersenne-Twister")
  defer_call(quote(suppressWarnings(RNGkind(.))), old_kind[[1L]])

  set.seed(8)
  invisible(fwb(d, lm_stat, R = 8, verbose = FALSE) |> futurize::futurize())

  expect_identical(RNGkind()[[1L]], "Mersenne-Twister")
})

test_that("`futurize()` replicates do not match an un-futurized call", {
  skip_if_no_futurize()

  d <- fwb_test_df()

  local_plan("multisession", workers = 2L)

  #Worth pinning explicitly, because it is the one place where `futurize()` is *not*
  #a transparent wrapper. `RNGkind("L'Ecuyer-CMRG")` re-initializes `.Random.seed`
  #even when L'Ecuyer-CMRG is already in force, so the stream *futurize()* hands to
  #`fwb()` is a deterministic function of the seed but not the seed itself. That
  #holds even with `simple = FALSE`, where every weight is drawn in the parent.
  #
  #The consequence for users: adding or removing `|> futurize()` changes the
  #replicates. Reproducing a published result means reproducing whether the pipe
  #was there.
  set.seed(4, "L'Ecuyer-CMRG")
  ref <- fwb(d, lm_stat, R = 20, verbose = FALSE, simple = FALSE)

  set.seed(4, "L'Ecuyer-CMRG")
  fz <- fwb(d, lm_stat, R = 20, verbose = FALSE, simple = FALSE) |>
    futurize::futurize()

  expect_different_t(fz, ref)

  #But `futurize()` is still internally consistent: with the weights drawn up
  #front, its replicates are invariant to the plan and the worker count.
  set.seed(4, "L'Ecuyer-CMRG")
  fz2 <- fwb(d, lm_stat, R = 20, verbose = FALSE, simple = FALSE) |>
    futurize::futurize()

  expect_same_t(fz2, fz)

  local_plan("multisession", workers = max_test_workers())

  set.seed(4, "L'Ecuyer-CMRG")
  fz3 <- fwb(d, lm_stat, R = 20, verbose = FALSE, simple = FALSE) |>
    futurize::futurize()

  expect_same_t(fz3, fz)
})

test_that("`futurize()` works with every weight type", {
  skip_if_no_futurize()

  d <- fwb_test_df()

  local_plan("multisession", workers = 2L)

  for (wt in c("exp", "multinom", "poisson", "mammen", "beta", "power")) {
    set.seed(12, "L'Ecuyer-CMRG")
    a <- fwb(d, lm_stat, R = 12, verbose = FALSE, wtype = wt) |>
      futurize::futurize()

    set.seed(12, "L'Ecuyer-CMRG")
    b <- fwb(d, lm_stat, R = 12, verbose = FALSE, wtype = wt) |>
      futurize::futurize()

    expect_identical(a[["wtype"]], wt)
    expect_same_t(b, a, info = wt)
  }
})

test_that("`futurize()` works with clusters and strata", {
  skip_if_no_futurize()

  d <- fwb_test_df()

  local_plan("multisession", workers = 2L)

  set.seed(13, "L'Ecuyer-CMRG")
  expect_no_error({
    a <- fwb(d, lm_stat, R = 16, verbose = FALSE, cluster = g) |>
      futurize::futurize()
  })

  set.seed(13, "L'Ecuyer-CMRG")
  b <- fwb(d, lm_stat, R = 16, verbose = FALSE, cluster = g) |>
    futurize::futurize()

  expect_same_t(b, a)
  expect_length(a[["cluster"]], nrow(d))

  #`simple = FALSE` because stratified sampling with `simple = TRUE` is broken for
  #the continuous weight types; see `test-strata.R`.
  set.seed(13, "L'Ecuyer-CMRG")
  expect_no_error({
    sa <- fwb(d, lm_stat, R = 16, verbose = FALSE, strata = s, simple = FALSE) |>
      futurize::futurize()
  })

  set.seed(13, "L'Ecuyer-CMRG")
  sb <- fwb(d, lm_stat, R = 16, verbose = FALSE, strata = s, simple = FALSE) |>
    futurize::futurize()

  expect_same_t(sb, sa)
})

test_that("`futurize()` gives reproducible `vcovFWB()` results", {
  skip_if_no_futurize()

  d <- fwb_test_df()
  fit <- lm(y ~ x1 + x2, data = d)

  local_plan("multisession", workers = 2L)

  set.seed(14, "L'Ecuyer-CMRG")
  expect_no_error({
    a <- vcovFWB(fit, R = 20) |> futurize::futurize()
  })

  set.seed(14, "L'Ecuyer-CMRG")
  b <- vcovFWB(fit, R = 20) |> futurize::futurize()

  expect_equal(a, b, tolerance = fwb_eps())
  expect_identical(dimnames(a), list(names(coef(fit)), names(coef(fit))))

  local_plan("sequential")

  set.seed(14, "L'Ecuyer-CMRG")
  s <- vcovFWB(fit, R = 20) |> futurize::futurize()

  expect_equal(s, a, tolerance = fwb_eps())
})

test_that("`futurize()` output supports the usual downstream methods", {
  skip_if_no_futurize()

  d <- fwb_test_df()

  local_plan("multisession", workers = 2L)

  #`futurize()` rewrites the call, so anything that reads `$call` -- `print()`,
  #and `fwb.array()` through the `cl` attribute -- sees a call the user never
  #wrote. These must all still work.
  set.seed(15, "L'Ecuyer-CMRG")
  f <- fwb(d, lm_stat, R = 80, verbose = FALSE) |> futurize::futurize()

  expect_output(print(f), "FRACTIONAL WEIGHTED BOOTSTRAP")

  expect_no_error({
    s <- summary(f)
  })

  expect_s3_class(s, "summary.fwb")

  expect_no_error({
    a <- fwb.array(f)
  })

  expect_identical(dim(a), c(80L, nrow(d)))
})
