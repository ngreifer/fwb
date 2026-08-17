#Replicability of `fwb()` and `vcovFWB()` across parallel backends.
#
#Since version 0.7.0 the contract is a single sentence: **the replicates depend on the
#seed and nothing else.** Neither `cl`, nor the worker count, nor `verbose`, nor the
#caller's RNG kind, nor randomness inside `statistic` can change them.
#
#That holds because the weights are never drawn from whatever stream a backend happens
#to hand a worker. `simple = FALSE` draws the whole matrix in the calling session before
#anything is dispatched; `simple = TRUE` records one L'Ecuyer stream per replicate up
#front and each replicate seeds its own weight draw from its own stream. Either way the
#backend contributes no randomness.
#
#These tests are the guard on that claim, and they are deliberately repetitive about it:
#a regression here would be invisible in ordinary use and would surface only as results
#that could not be reproduced on another machine.

#Every configuration worth running `fwb()` under, as a list of argument lists. `setup`
#runs first and may `skip()`; it returns extra arguments to splice into the call.
backend_configs <- function() {
  list(
    list(label = "sequential",
         setup = function() list()),
    list(label = "verbose",
         setup = function() list(verbose = TRUE)),
    list(label = "cl = 2",
         setup = function() {
           skip_if_no_forking()
           list(cl = 2)
         }),
    list(label = sprintf("cl = %d", max_test_workers()),
         setup = function() {
           skip_if_no_forking()
           list(cl = max_test_workers())
         }),
    list(label = "future, 2 workers",
         setup = function() {
           local_plan("multisession", workers = 2L, .env = parent.frame(2L))
           list(cl = "future")
         }),
    list(label = sprintf("future, %d workers", max_test_workers()),
         setup = function() {
           local_plan("multisession", workers = max_test_workers(),
                      .env = parent.frame(2L))
           list(cl = "future")
         }),
    list(label = "future, verbose",
         setup = function() {
           local_plan("multisession", workers = 2L, .env = parent.frame(2L))
           list(cl = "future", verbose = TRUE)
         }),
    list(label = "cluster object",
         setup = function() {
           list(cl = local_cluster(2L, .env = parent.frame(2L)))
         })
  )
}

#Run `fwb()` under one configuration from the same seed, discarding printed output so a
#progress bar does not litter the reporter.
run_under <- function(config, ..., seed = 1234) {
  extra <- config[["setup"]]()

  args <- utils::modifyList(list(verbose = FALSE), extra)

  set.seed(seed, "L'Ecuyer-CMRG")

  out <- NULL
  utils::capture.output({
    out <- do.call(fwb, c(list(fwb_test_df(), lm_stat, R = 24L), args, list(...)))
  })

  out
}

test_that("the replicates are the same under every backend", {
  skip_on_cran()

  for (simple in c(TRUE, FALSE)) {
    configs <- backend_configs()

    ref <- run_under(configs[[1L]], simple = simple)

    for (config in configs[-1L]) {
      out <- run_under(config, simple = simple)

      expect_same_t(out, ref,
                    info = sprintf("simple = %s, %s", simple, config[["label"]]))
    }
  }
})

test_that("the replicates are the same under every backend, for every weight type", {
  skip_on_cran()

  #The weight types differ in how they consume the stream -- `"multinom"` in particular
  #fills its batch column-major -- so backend invariance has to be checked per type
  #rather than assumed from the default.
  local_plan("multisession", workers = 2L)

  cl <- local_cluster(2L)

  d <- fwb_test_df()

  for (wt in c("exp", "multinom", "poisson", "mammen", "beta", "power")) {
    for (simple in c(TRUE, FALSE)) {
      set.seed(99, "L'Ecuyer-CMRG")
      ref <- fwb(d, lm_stat, R = 16L, verbose = FALSE, wtype = wt, simple = simple)

      set.seed(99, "L'Ecuyer-CMRG")
      fut <- fwb(d, lm_stat, R = 16L, verbose = FALSE, wtype = wt, simple = simple,
                 cl = "future")

      set.seed(99, "L'Ecuyer-CMRG")
      clus <- fwb(d, lm_stat, R = 16L, verbose = FALSE, wtype = wt, simple = simple,
                  cl = cl)

      info <- sprintf("%s, simple = %s", wt, simple)

      expect_same_t(fut, ref, info = info)
      expect_same_t(clus, ref, info = info)
    }
  }
})

test_that("the replicates are the same under every backend, with clusters and strata", {
  skip_on_cran()

  local_plan("multisession", workers = 2L)

  d <- fwb_test_df()

  for (extra in list(list(cluster = quote(g)),
                     list(strata = quote(s)),
                     list(cluster = quote(g), strata = quote(s)))) {
    for (simple in c(TRUE, FALSE)) {
      call_args <- c(list(d, lm_stat, R = 16L, verbose = FALSE, simple = simple),
                     extra)

      set.seed(31, "L'Ecuyer-CMRG")
      ref <- do.call(fwb, call_args)

      set.seed(31, "L'Ecuyer-CMRG")
      fut <- do.call(fwb, c(call_args, list(cl = "future")))

      expect_same_t(fut, ref,
                    info = sprintf("%s, simple = %s",
                                   paste(names(extra), collapse = "+"), simple))
    }
  }
})

test_that("a plain `set.seed()` is enough, under any RNG kind", {
  #Before 0.7.0 the weights were drawn inside the workers, so reproducing a parallel run
  #meant `set.seed(n, "L'Ecuyer-CMRG")` -- and with a `cluster` object, not `set.seed()`
  #at all but `parallel::clusterSetRNGStream()`. Neither is needed now: the streams are
  #derived in the calling session, and their derivation does not depend on the caller's
  #kind.
  skip_on_cran()

  d <- fwb_test_df()

  old_kind <- RNGkind("Mersenne-Twister")
  defer_call(quote(suppressWarnings(RNGkind(.))), old_kind[[1L]])

  for (cl_arg in list(NULL, 2)) {
    if (!is.null(cl_arg)) {
      skip_if_no_forking()
    }

    RNGkind("Mersenne-Twister")
    set.seed(6)
    a <- fwb(d, lm_stat, R = 20L, verbose = FALSE, cl = cl_arg)

    RNGkind("Mersenne-Twister")
    set.seed(6)
    b <- fwb(d, lm_stat, R = 20L, verbose = FALSE, cl = cl_arg)

    expect_same_t(b, a, info = format(cl_arg))

    RNGkind("Mersenne-Twister")
    set.seed(7)
    other <- fwb(d, lm_stat, R = 20L, verbose = FALSE, cl = cl_arg)

    expect_different_t(other, a, info = format(cl_arg))
  }

  #The caller's RNG kind is left as it was, even though the recorded streams are
  #L'Ecuyer and installing one would otherwise switch it.
  expect_identical(RNGkind()[[1L]], "Mersenne-Twister")

  #And the same seed under a different kind gives a different answer, which is the
  #price of not requiring one: the streams are derived from the caller's state, and
  #that state is kind-specific.
  RNGkind("Mersenne-Twister")
  set.seed(6)
  mt <- fwb(d, lm_stat, R = 20L, verbose = FALSE)

  set.seed(6, "L'Ecuyer-CMRG")
  lec <- fwb(d, lm_stat, R = 20L, verbose = FALSE)

  expect_different_t(lec, mt)
})

test_that("successive calls without a new seed give different replicates", {
  d <- fwb_test_df()

  #Deriving the streams consumes exactly one draw from the caller's stream, so a second
  #call proceeds from a different state. Without that, `fwb()` would silently return the
  #same bootstrap twice.
  set.seed(11, "L'Ecuyer-CMRG")
  a <- fwb(d, lm_stat, R = 16L, verbose = FALSE)
  b <- fwb(d, lm_stat, R = 16L, verbose = FALSE)

  expect_different_t(b, a)

  #And the pair is itself reproducible.
  set.seed(11, "L'Ecuyer-CMRG")
  a2 <- fwb(d, lm_stat, R = 16L, verbose = FALSE)
  b2 <- fwb(d, lm_stat, R = 16L, verbose = FALSE)

  expect_same_t(a2, a)
  expect_same_t(b2, b)
})

test_that("`parallel::clusterSetRNGStream()` no longer affects the result", {
  skip_on_cran()

  d <- fwb_test_df()

  cl <- local_cluster(2L)

  #The workers' streams are never used to draw weights, so setting them does nothing --
  #which is why `vignette("fwb-rep")`'s old Case 5 advice no longer applies. Two runs
  #seeded only by `clusterSetRNGStream()` differ, and two runs seeded by `set.seed()`
  #agree.
  parallel::clusterSetRNGStream(cl, 4321)
  a <- fwb(d, lm_stat, R = 24L, verbose = FALSE, cl = cl)

  parallel::clusterSetRNGStream(cl, 4321)
  b <- fwb(d, lm_stat, R = 24L, verbose = FALSE, cl = cl)

  expect_different_t(b, a)

  set.seed(11, "L'Ecuyer-CMRG")
  s1 <- fwb(d, lm_stat, R = 24L, verbose = FALSE, cl = cl)

  set.seed(11, "L'Ecuyer-CMRG")
  s2 <- fwb(d, lm_stat, R = 24L, verbose = FALSE, cl = cl)

  expect_same_t(s2, s1)
})

test_that("randomness in `statistic` does not disturb the weights", {
  skip_on_cran()

  d <- fwb_test_df()

  #Each replicate's weight draw is seeded from its own recorded stream, so it cannot be
  #shifted by anything `statistic` drew in an earlier replicate -- the failure mode that
  #made the old design give up on this case entirely.
  set.seed(21, "L'Ecuyer-CMRG")
  a <- fwb(d, random_lm_stat, R = 20L, verbose = FALSE)

  set.seed(21, "L'Ecuyer-CMRG")
  b <- fwb(d, random_lm_stat, R = 20L, verbose = FALSE)

  expect_same_t(a, b)

  #Including in parallel, and with the same weights as the sequential run: the weights
  #match a run whose `statistic` drew nothing at all.
  local_plan("multisession", workers = 2L)

  set.seed(21, "L'Ecuyer-CMRG")
  par <- fwb(d, random_lm_stat, R = 20L, verbose = FALSE, cl = "future")

  expect_same_t(par, a)

  set.seed(21, "L'Ecuyer-CMRG")
  w_random <- fwb(d, function(data, w) as.numeric(w) + 0 * rnorm(1L),
                  R = 20L, verbose = FALSE)

  set.seed(21, "L'Ecuyer-CMRG")
  w_plain <- fwb(d, w_identity_stat, R = 20L, verbose = FALSE)

  expect_same_t(w_random, w_plain)
})

test_that("`simple` changes the replicates but not their distribution", {
  d <- fwb_test_df()

  #`simple = FALSE` draws the whole matrix in one call, which is what keeps
  #`wtype = "multinom"` identical to `boot::boot(., stype = "f")`; `simple = TRUE` draws
  #each replicate from its own stream. The two consume the stream differently, so they
  #no longer agree -- worth pinning, because before 0.7.0 they did for the continuous
  #weight types and someone may remember that.
  set.seed(41, "L'Ecuyer-CMRG")
  a <- fwb(d, w_identity_stat, R = 400L, verbose = FALSE, simple = TRUE)

  set.seed(41, "L'Ecuyer-CMRG")
  b <- fwb(d, w_identity_stat, R = 400L, verbose = FALSE, simple = FALSE)

  expect_different_t(a, b)

  #Both are draws from the same weight distribution, though.
  expect_equal(mean(a[["t"]]), 1, tolerance = 0.01)
  expect_equal(mean(b[["t"]]), 1, tolerance = 0.01)
  expect_equal(var(as.vector(a[["t"]])), var(as.vector(b[["t"]])),
               tolerance = 0.1)
})

test_that("`verbose` defaults to `FALSE` only when parallelizing", {
  d <- fwb_test_df()

  #`verbose = NULL` resolves to `guess_num_workers(cl) == 1`, so the sequential default
  #still prints a progress bar and the parallel default does not.
  set.seed(1, "L'Ecuyer-CMRG")
  out <- utils::capture.output({
    f <- fwb(d, lm_stat, R = 8L)
  }, type = "output")

  expect_true(any(nzchar(out)))

  skip_on_cran()
  local_plan("multisession", workers = 2L)

  set.seed(1, "L'Ecuyer-CMRG")
  out <- utils::capture.output({
    f <- fwb(d, lm_stat, R = 8L, cl = "future")
  }, type = "output")

  expect_false(any(nzchar(out)))
})

test_that("`fwb()` no longer warns about the RNG kind", {
  skip_on_cran()
  skip_if_no_forking()

  d <- fwb_test_df()

  #The warning told users to `set.seed(n, "L'Ecuyer-CMRG")` before parallelizing, which
  #was necessary while the workers drew the weights. It is not any more, and telling
  #people to do something that no longer matters is worse than saying nothing.
  old_kind <- RNGkind("Mersenne-Twister")
  defer_call(quote(suppressWarnings(RNGkind(.))), old_kind[[1L]])

  set.seed(1)
  expect_no_warning({
    f <- fwb(d, lm_stat, R = 8L, verbose = FALSE, cl = 2)
  }, message = "not suitable for parallelization")

  set.seed(1)
  expect_no_warning({
    v <- vcovFWB(lm(y ~ x1 + x2, data = d), R = 8L, cl = 2)
  }, message = "not suitable for parallelization")
})

test_that("what each `simple` setting records is enough to recover the weights", {
  d <- fwb_test_df()

  #Both settings record what they need in the `"seeds"` attribute: `simple = TRUE` one
  #stream per replicate, `simple = FALSE` the single state its one batch call was made
  #from. `fwb.array()` reads whichever is there -- `test-fwb-array.R` covers the round
  #trip, this covers the bookkeeping.
  set.seed(31, "L'Ecuyer-CMRG")
  fs <- fwb(d, lm_stat, R = 12L, verbose = FALSE, simple = TRUE)

  seeds <- attr(fs, "seeds", TRUE)

  expect_identical(dim(seeds), c(12L, 7L))
  expect_identical(anyDuplicated(seeds), 0L)

  #Each row is an L'Ecuyer-CMRG stream, and consecutive rows are consecutive streams.
  expect_identical(seeds[2L, ], parallel::nextRNGStream(seeds[1L, ]))

  set.seed(31, "L'Ecuyer-CMRG")
  fb <- fwb(d, lm_stat, R = 12L, verbose = FALSE, simple = FALSE)

  seed <- attr(fb, "seeds", TRUE)

  expect_type(seed, "integer")
  expect_false(is.matrix(seed))

  #And that one state is enough: restoring it and re-drawing in a single call reproduces
  #the weights the bootstrap used.
  gen <- fwb:::make_gen_weights(fb[["wtype"]])

  w <- fwb:::with_seed_preserved(gen(nrow(d), fb[["R"]], NULL), new_seed = seed)

  set.seed(31, "L'Ecuyer-CMRG")
  fb_w <- fwb(d, w_identity_stat, R = 12L, verbose = FALSE, simple = FALSE)

  expect_equal(unname(fb_w[["t"]]), unname(w), tolerance = fwb_eps())
})

test_that("`with_seed_preserved()` leaves the stream untouched", {
  set.seed(31, "L'Ecuyer-CMRG")
  before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)

  invisible(fwb:::with_seed_preserved(runif(5L)))

  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   before)

  #Including when a foreign stream is installed, which also switches the RNG kind.
  set.seed(31, "L'Ecuyer-CMRG")
  lecuyer <- get(".Random.seed", envir = globalenv(), inherits = FALSE)

  old_kind <- RNGkind("Mersenne-Twister")
  defer_call(quote(suppressWarnings(RNGkind(.))), old_kind[[1L]])

  set.seed(3)
  before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)

  invisible(fwb:::with_seed_preserved(runif(5L), new_seed = lecuyer))

  expect_identical(RNGkind()[[1L]], "Mersenne-Twister")
  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   before)
})

test_that("`vcovFWB()` is the same under every backend", {
  skip_on_cran()

  d <- fwb_test_df()
  fit <- lm(y ~ x1 + x2, data = d)

  set.seed(1234, "L'Ecuyer-CMRG")
  ref <- vcovFWB(fit, R = 24L)

  #`vcovFWB()` draws its weights the same way `fwb(simple = TRUE)` does, so it inherits
  #the same guarantee.
  if (.Platform$OS.type != "windows") {
    for (k in unique(c(2L, max_test_workers()))) {
      set.seed(1234, "L'Ecuyer-CMRG")
      expect_equal(vcovFWB(fit, R = 24L, cl = k), ref, tolerance = fwb_eps(),
                   info = paste("cl =", k))
    }
  }

  for (workers in unique(c(2L, max_test_workers()))) {
    local_plan("multisession", workers = workers)

    set.seed(1234, "L'Ecuyer-CMRG")
    expect_equal(vcovFWB(fit, R = 24L, cl = "future"), ref,
                 tolerance = fwb_eps(), info = paste("future,", workers))
  }

  cl <- local_cluster(2L)

  set.seed(1234, "L'Ecuyer-CMRG")
  expect_equal(vcovFWB(fit, R = 24L, cl = cl), ref, tolerance = fwb_eps())

  set.seed(1234, "L'Ecuyer-CMRG")
  utils::capture.output({
    verb <- vcovFWB(fit, R = 24L, verbose = TRUE)
  })

  expect_equal(verb, ref, tolerance = fwb_eps())
})

test_that("`vcovFWB()` agrees with `fwb(simple = TRUE)` on every weight type", {
  d <- fwb_test_df()

  #The documented relationship between the two functions: `vcovFWB()` is `fwb()` plus
  #`vcov()` on the same weights. Both now seed per replicate from the same derivation,
  #which is what keeps that true.
  specs <- list(
    list(fit = lm(y ~ x1 + x2, data = d), stat = lm_stat),
    list(fit = glm(yb ~ x1 + x2, data = d, family = quasibinomial()),
         stat = glm_stat))

  for (spec in specs) {
    for (wt in c("exp", "multinom", "poisson", "mammen", "beta", "power")) {
      set.seed(101, "L'Ecuyer-CMRG")
      v <- vcovFWB(spec[["fit"]], R = 40L, wtype = wt)

      set.seed(101, "L'Ecuyer-CMRG")
      f <- fwb(d, spec[["stat"]], R = 40L, verbose = FALSE, wtype = wt,
               simple = TRUE)

      expect_equal(v, vcov(f), tolerance = 1e-6,
                   info = paste(class(spec[["fit"]])[1L], wt))
    }
  }
})

test_that("`vcovFWB()` draws each cluster dimension independently", {
  d <- fwb_test_df()
  fit <- lm(y ~ x1 + x2, data = d)

  #Multi-way clustering runs the replicates once per dimension. Each dimension takes a
  #fresh set of recorded streams, so the dimensions are independent -- reusing one set
  #would correlate the terms that the inclusion-exclusion sum is supposed to combine.
  set.seed(6, "L'Ecuyer-CMRG")
  one_way <- vcovFWB(fit, cluster = ~ g, R = 30L)

  set.seed(6, "L'Ecuyer-CMRG")
  two_way <- vcovFWB(fit, cluster = ~ g + x2, R = 30L)

  expect_not_equal(one_way, two_way)

  set.seed(6, "L'Ecuyer-CMRG")
  two_way_again <- vcovFWB(fit, cluster = ~ g + x2, R = 30L)

  expect_equal(two_way_again, two_way, tolerance = fwb_eps())
})

test_that("`vcovFWB()` works when the RNG has not been used yet", {
  #`vcovFWB()` reaches `with_seed_preserved()` before anything has drawn from the RNG,
  #and `.Random.seed` does not exist in a session that has never used it -- as in
  #`Rscript --vanilla`, or a fresh `callr`/`targets` worker. The guard now lives in
  #`with_seed_preserved()` itself, so both `fwb()` and `vcovFWB()` have it.
  #
  #Deleting `.Random.seed` reproduces that state in-process. It is restored first thing
  #on exit so a failure here cannot leave the rest of the run without a stream.
  old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  defer_call(quote(assign(".Random.seed", ., envir = globalenv())), old_seed)

  rm(".Random.seed", envir = globalenv())

  d <- fwb_test_df()

  expect_no_error({
    v <- vcovFWB(lm(y ~ x1 + x2, data = d), R = 10L)
  })

  expect_identical(dim(v), c(3L, 3L))

  rm(".Random.seed", envir = globalenv())

  expect_no_error({
    f <- fwb(d, lm_stat, R = 10L, verbose = FALSE)
  })

  expect_identical(dim(f[["t"]]), c(10L, 3L))
})
