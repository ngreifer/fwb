#`fwb.array()` re-generates the bootstrap weights from the saved `seed`. Every test
#here is the same round trip: run `fwb()` with a `statistic` that *returns its own
#weights*, then check that `fwb.array()` reproduces them exactly.
#
#That comparison is worth having as its own file because `fwb.array()` is not only a
#convenience -- `compute_ci()` calls it to build the empirical influence values for
#BCa intervals. A silent mismatch there does not error; it produces a plausible-
#looking interval computed from the wrong influence function.

expect_array_round_trip <- function(..., R = 24L, n = 60L, tolerance = fwb_eps()) {
  d <- fwb_test_df(n)

  set.seed(64, "L'Ecuyer-CMRG")
  f <- fwb(d, w_identity_stat, R = R, verbose = FALSE, ...)

  a <- fwb.array(f)

  expect_identical(dim(a), c(R, n))
  expect_equal(unname(f[["t"]]), unname(a), tolerance = tolerance)
}

test_that("`fwb.array()` recovers the weights for the continuous weight types", {
  for (wt in c("exp", "mammen", "beta", "power")) {
    for (simple in c(TRUE, FALSE)) {
      expect_array_round_trip(wtype = wt, simple = simple)
    }
  }
})

test_that("`fwb.array()` recovers the weights for the integer weight types", {
  #`"poisson"` fills its matrix row-wise, so batch and per-replicate generation
  #consume the RNG identically and both values of `simple` round-trip.
  expect_array_round_trip(wtype = "poisson", simple = TRUE)
  expect_array_round_trip(wtype = "poisson", simple = FALSE)

  #`"multinom"` with its default `simple = FALSE`.
  expect_array_round_trip(wtype = "multinom", simple = FALSE)
})

test_that("`fwb.array()` recovers `\"multinom\"` weights when `simple = TRUE`", {
  #`fwb.array()` chooses between batch and per-replicate re-generation. The batch call
  #is only equivalent for generators whose `R`-row fill consumes the stream the way `R`
  #one-row calls do, which `"multinom"` does not: it fills column-major, because that
  #is what makes it reproduce `boot::boot(., stype = "f")`. So for `"multinom"` with
  #`simple = TRUE` the batch path returns a reshuffle of the weights actually used --
  #right dimensions, right marginals, wrong rows, and no error.
  expect_array_round_trip(wtype = "multinom", simple = TRUE)

  #Same with strata, which uses the other branch of the generator.
  d <- fwb_test_df()

  set.seed(64, "L'Ecuyer-CMRG")
  f <- fwb(d, w_identity_stat, R = 24L, verbose = FALSE, wtype = "multinom",
           simple = TRUE, strata = s)

  expect_equal(unname(f[["t"]]), unname(fwb.array(f)), tolerance = fwb_eps())
})

test_that("`fwb.array()` recovers the weights with a factor `strata`", {
  #`simple = TRUE` with strata is tested in `test-strata.R`, where it currently
  #errors; `simple = FALSE` is the working path.
  for (wt in c("exp", "mammen", "beta", "power", "multinom")) {
    d <- fwb_test_df()

    set.seed(65, "L'Ecuyer-CMRG")
    f <- fwb(d, w_identity_stat, R = 24L, verbose = FALSE, wtype = wt,
             simple = FALSE, strata = s)

    a <- fwb.array(f)

    expect_equal(unname(f[["t"]]), unname(a), tolerance = fwb_eps(), info = wt)
  }
})

test_that("`fwb.array()` recovers the weights with a non-factor `strata`", {
  #`fwb()` stores `strata` exactly as supplied and `fwb.array()` hands that stored
  #vector straight to the generator, which decides whether to stratify by calling
  #`nlevels()`. For a character or numeric vector `nlevels()` is 0, so the
  #re-generated weights are *unstratified* -- they no longer sum to the stratum
  #sizes, and they are not the weights the bootstrap used.
  #
  #`fwb()` itself is unaffected: it converts to a factor internally. Only the
  #recovery path is wrong, which means the BCa influence values are computed from
  #the wrong weights whenever `strata` was not already a factor.
  #
  #The `factor()` conversion has to happen above the early `return()` for the batch path,
  #which is the one taken whenever `cl` is `NULL` -- for either value of `simple`.
  d <- fwb_test_df()
  d[["s_chr"]] <- as.character(d[["s"]])

  for (simple in c(TRUE, FALSE)) {
    set.seed(65, "L'Ecuyer-CMRG")
    f <- fwb(d, w_identity_stat, R = 24L, verbose = FALSE, simple = simple,
             strata = s_chr)

    expect_equal(unname(f[["t"]]), unname(fwb.array(f)), tolerance = fwb_eps(),
                 info = paste("simple =", simple))
  }
})

test_that("`fwb.array()` recovers cluster weights", {
  #With clusters, `fwb()` draws one weight per cluster and spreads it over the cluster's
  #members, so a recovered row has at most `nlevels(cluster)` distinct values.
  d <- fwb_test_df()

  for (simple in c(TRUE, FALSE)) {
    set.seed(66, "L'Ecuyer-CMRG")
    f <- fwb(d, w_identity_stat, R = 24L, verbose = FALSE, cluster = g,
             simple = simple)

    a <- fwb.array(f)

    expect_identical(dim(a), c(24L, nrow(d)))
    expect_equal(unname(f[["t"]]), unname(a), tolerance = fwb_eps(),
                 info = paste("simple =", simple))

    #Members of a cluster share a weight.
    expect_true(all(apply(a, 1L, function(w) {
      all(tapply(w, d[["g"]], function(z) fwb:::all_the_same(z)))
    })))

    expect_true(all(apply(a, 1L, function(w) {
      length(unique(round(w, 10L))) <= nlevels(d[["g"]])
    })))
  }

  #And with strata as well, where `fwb()` collapses `strata` to one entry per cluster
  #before drawing -- a collapse `fwb.array()` has to reproduce or it would hand the
  #generator the wrong number of columns.
  for (simple in c(TRUE, FALSE)) {
    set.seed(66, "L'Ecuyer-CMRG")
    f <- fwb(d, w_identity_stat, R = 24L, verbose = FALSE, cluster = g,
             strata = s, simple = simple)

    expect_equal(unname(f[["t"]]), unname(fwb.array(f)), tolerance = fwb_eps(),
                 info = paste("simple =", simple))
  }
})

test_that("`fwb.array()` recovers the weights across parallel backends", {
  skip_on_cran()

  d <- fwb_test_df()

  #Recovery reads the recorded streams, so it does not matter how the original run was
  #parallelized -- and `fwb.array()` needs no backend of its own.
  set.seed(67, "L'Ecuyer-CMRG")
  ref <- fwb(d, w_identity_stat, R = 24L, verbose = FALSE)

  expect_equal(unname(ref[["t"]]), unname(fwb.array(ref)), tolerance = fwb_eps())

  if (.Platform$OS.type != "windows") {
    set.seed(67, "L'Ecuyer-CMRG")
    f <- fwb(d, w_identity_stat, R = 24L, verbose = FALSE, cl = 2)

    expect_equal(unname(f[["t"]]), unname(fwb.array(f)), tolerance = fwb_eps())
    expect_same_t(f, ref)
  }

  local_plan("multisession", workers = 2L)

  for (simple in c(TRUE, FALSE)) {
    set.seed(67, "L'Ecuyer-CMRG")
    f <- fwb(d, w_identity_stat, R = 24L, verbose = FALSE, cl = "future",
             simple = simple)

    expect_equal(unname(f[["t"]]), unname(fwb.array(f)), tolerance = fwb_eps(),
                 info = paste("simple =", simple))
  }
})

test_that("`fwb.array()` needs no RNG setup from the caller", {
  skip_on_cran()

  d <- fwb_test_df()

  cl <- local_cluster(2L)

  #`vignette("fwb-rep")` used to tell users to call `clusterSetRNGStream()` with the
  #original seed before `fwb.ci()`/`summary()`, because recovery re-ran the map on the
  #workers. Recovery now reads the recorded streams in the calling session, so repeated
  #calls agree and nothing has to be set up first.
  set.seed(2468, "L'Ecuyer-CMRG")
  f <- fwb(d, w_identity_stat, R = 24L, verbose = FALSE, cl = cl)

  expect_equal(unname(f[["t"]]), unname(fwb.array(f)), tolerance = fwb_eps())
  expect_identical(fwb.array(f), fwb.array(f))

  #Even after the caller's stream has moved on.
  invisible(runif(10L))

  expect_equal(unname(f[["t"]]), unname(fwb.array(f)), tolerance = fwb_eps())

  #And after the workers are gone, because the object never held on to them.
  parallel::stopCluster(cl)

  expect_equal(unname(f[["t"]]), unname(fwb.array(f)), tolerance = fwb_eps())

  #Likewise after a round trip through `saveRDS()`, which a live `cluster` object in the
  #result used to make impossible.
  path <- tempfile(fileext = ".rds")
  saveRDS(f, path)
  on.exit(unlink(path), add = TRUE)

  f2 <- readRDS(path)

  expect_equal(unname(f2[["t"]]), unname(fwb.array(f2)), tolerance = fwb_eps())
})

test_that("`fwb.array()` leaves the caller's RNG state alone", {
  d <- fwb_test_df()

  set.seed(68, "L'Ecuyer-CMRG")
  f <- fwb(d, lm_stat, R = 16L, verbose = FALSE)

  before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)

  invisible(fwb.array(f))

  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   before)

  #Also true when the array is re-generated through a backend.
  skip_on_cran()
  local_plan("multisession", workers = 2L)

  set.seed(68, "L'Ecuyer-CMRG")
  f <- fwb(d, lm_stat, R = 16L, verbose = FALSE, cl = "future")

  before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)

  invisible(fwb.array(f))

  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   before)
})

test_that("`fwb.array()` recovers the weights even when `statistic` is random", {
  skip_on_cran()

  d <- fwb_test_df()

  #This is the case the replay-the-stream design could not handle: with `simple = TRUE`
  #and a `statistic` that draws its own random numbers, the weight draws and the
  #statistic's draws were interleaved in an order that could not be replayed, so
  #`fwb.array()` warned and returned the wrong weights, and BCa intervals were refused.
  #
  #Seeding each replicate's weights from its own recorded stream removes the
  #interleaving: the weights are a function of the replicate index and nothing else.
  random_w_stat <- function(data, w) as.numeric(w) + 0 * rnorm(1L)

  set.seed(69, "L'Ecuyer-CMRG")
  f <- fwb(d, random_w_stat, R = 16L, verbose = FALSE)

  expect_no_warning({
    a <- fwb.array(f)
  })

  expect_equal(unname(f[["t"]]), unname(a), tolerance = fwb_eps())

  if (.Platform$OS.type != "windows") {
    set.seed(69, "L'Ecuyer-CMRG")
    fp <- fwb(d, random_w_stat, R = 16L, verbose = FALSE, cl = 2)

    expect_no_warning({
      ap <- fwb.array(fp)
    })

    expect_equal(unname(fp[["t"]]), unname(ap), tolerance = fwb_eps())
  }

  #So BCa intervals are available here too, which they previously were not.
  d40 <- fwb_test_df(40L)

  set.seed(69, "L'Ecuyer-CMRG")
  fr <- fwb(d40, function(data, w) lm_stat(data, w) + rnorm(3L, sd = 1e-7),
            R = 200L, verbose = FALSE)

  expect_no_error({
    ci <- fwb.ci(fr, type = "bca", index = 2L)
  })

  expect_s3_class(ci, "fwbci")
})

test_that("`fwb.array()` rejects objects it cannot handle", {
  d <- fwb_test_df()

  expect_err(fwb.array(), "must be supplied")
  expect_err(fwb.array(d), "must inherit from class")

  set.seed(70, "L'Ecuyer-CMRG")
  f <- fwb(d, lm_stat, R = 8L, verbose = FALSE)

  #`fwb.array()` also accepts `<boot>` objects, delegating to `boot::boot.array()`.
  skip_if_not_installed("boot")

  set.seed(70)
  b <- boot::boot(d, function(dat, i) coef(lm(y ~ x1 + x2, data = dat[i, ])),
                  R = 8L)

  expect_identical(fwb.array(b), boot::boot.array(b))
})

test_that("BCa intervals use the re-generated weights", {
  d <- fwb_test_df(40L)

  #BCa needs `R > n`, and it is the only consumer of `fwb.array()` inside the
  #package. Computing the acceleration by hand from the array and reproducing the
  #interval is what ties the two together: if `fwb.array()` started returning
  #different weights, this would break even though `fwb.ci()` itself did not change.
  set.seed(71, "L'Ecuyer-CMRG")
  f <- fwb(d, lm_stat, R = 200L, verbose = FALSE, simple = FALSE)

  w <- fwb.array(f)

  index <- 2L
  t <- f[["t"]][, index, drop = FALSE]

  empinf <- fwb:::empinf.reg(f, t = t)

  a_hand <- sum(empinf^3) / (6 * sum(empinf^2)^1.5)

  expect_equal(unname(a_hand),
               unname(sum(empinf^3) / (6 * sum(empinf^2)^1.5)),
               tolerance = fwb_eps())

  #The influence values must be centered and have one entry per unit.
  expect_identical(dim(empinf), c(nrow(d), 1L))
  expect_equal(mean(empinf), 0, tolerance = fwb_eps())

  #And they must be a function of the recovered weights, not of anything else.
  expinf2 <- as.matrix(.lm.fit(x = w / nrow(d), y = t)$coefficients)
  expinf2[] <- sweep(expinf2, 2L, colMeans(expinf2))

  expect_equal(unname(empinf), unname(expinf2), tolerance = fwb_eps())

  expect_no_error({
    ci <- fwb.ci(f, type = "bca", index = index)
  })

  expect_s3_class(ci, "fwbci")
})
