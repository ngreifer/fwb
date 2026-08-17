#Argument validation for `fwb()` and `vcovFWB()`.
#
#Bad input to a bootstrap is expensive to diagnose after the fact: `R` model fits
#have already run, and a mistake in `statistic` shows up as a wall of replicate
#estimates rather than as an error. Each expectation below pins one guard to the
#message it produces, so a guard that stops firing -- or starts firing on valid
#input -- fails here rather than in someone's analysis.

test_that("`fwb()` requires usable `data` and `statistic`", {
  d <- fwb_test_df()

  expect_err(fwb(), "must be supplied")
  expect_err(fwb(d), "must be supplied")
  expect_err(fwb(1:10, lm_stat, R = 5L), "must be a data frame")
  expect_err(fwb(d, "lm_stat", R = 5L), "must be a function")

  #A zero-row data frame has nothing to draw weights for.
  expect_err(fwb(d[0L, ], lm_stat, R = 5L), "must be present")
})

test_that("`fwb()` requires a positive whole-number `R`", {
  d <- fwb_test_df()

  expect_err(fwb(d, lm_stat, R = 0L, verbose = FALSE), "must be positive")
  expect_err(fwb(d, lm_stat, R = -5L, verbose = FALSE), "must be")
  expect_err(fwb(d, lm_stat, R = 2.5, verbose = FALSE), "must be")
  expect_err(fwb(d, lm_stat, R = c(5L, 10L), verbose = FALSE), "must be")

  expect_no_error({
    f <- fwb(d, lm_stat, R = 1L, verbose = FALSE)
  })

  expect_identical(f[["R"]], 1L)
  expect_identical(nrow(f[["t"]]), 1L)
})

test_that("`fwb()` rejects a `statistic` that does not return a numeric vector", {
  d <- fwb_test_df()

  expect_err(fwb(d, function(data, w) matrix(1:4, 2L), R = 5L, verbose = FALSE),
             "must be a numeric vector")
  expect_err(fwb(d, function(data, w) "a", R = 5L, verbose = FALSE),
             "must be a numeric vector")
  expect_err(fwb(d, function(data, w) list(1, 2), R = 5L, verbose = FALSE),
             "must be a numeric vector")

  #`NA` in the original sample means there is nothing to bootstrap around.
  expect_err(fwb(d, function(data, w) c(a = NA_real_), R = 5L, verbose = FALSE),
             "returned as NA in the original sample")

  #An error inside `statistic` is reported as such, with the original error attached
  #as the parent rather than swallowed.
  expect_err(fwb(d, function(data, w) stop("kaboom"), R = 5L, verbose = FALSE),
             "there was a problem running the function supplied to")
  expect_err(fwb(d, function(data, w) stop("kaboom"), R = 5L, verbose = FALSE),
             "kaboom")
})

test_that("`fwb()` warns about `NA` replicates without discarding them", {
  d <- fwb_test_df()

  #An occasional `NA` replicate is not fatal -- it may be a model that failed to
  #converge for one set of weights -- but it does break the intervals downstream, so
  #the user has to hear about it once.
  flaky <- function(data, w) {
    est <- coef(lm(y ~ x1, data = data, weights = w))

    if (w[[1L]] > 2) {
      is.na(est) <- TRUE
    }

    est
  }

  set.seed(93, "L'Ecuyer-CMRG")
  expect_wrn(f <- fwb(d, flaky, R = 40L, verbose = FALSE),
             "returned as NA")

  expect_true(anyNA(f[["t"]]))
  expect_identical(nrow(f[["t"]]), 40L)

  #`print()` still works, reporting the non-`NA` statistics.
  expect_output(print(f), "Bootstrap Statistics")
})

test_that("`fwb()` reports when every replicate of a statistic is `NA`", {
  d <- fwb_test_df()

  all_na <- function(data, w) {
    c(good = coef(lm(y ~ x1, data = data, weights = w))[[2L]],
      bad = if (all(w == 1)) 0 else NA_real_)
  }

  set.seed(94, "L'Ecuyer-CMRG")
  expect_wrn(f <- fwb(d, all_na, R = 20L, verbose = FALSE), "returned as NA")

  out <- collapse_ws(utils::capture.output(print(f)))

  expect_match(out, "WARNING: All values of bad* are NA", fixed = TRUE)
  expect_match(out, "good", fixed = TRUE)
})

test_that("`fwb()` validates `cluster` and `strata` lengths", {
  d <- fwb_test_df()

  expect_err(fwb(d, lm_stat, R = 5L, verbose = FALSE, cluster = 1:5),
             "must have length 60")
  expect_err(fwb(d, lm_stat, R = 5L, verbose = FALSE, strata = 1:5,
                 simple = FALSE),
             "must have length 60")

  #A matrix is not a valid membership vector.
  expect_err(fwb(d, lm_stat, R = 5L, verbose = FALSE,
                 cluster = matrix(1L, nrow(d), 2L)),
             "must be")
})

test_that("`fwb()` validates the remaining flags", {
  d <- fwb_test_df()

  expect_err(fwb(d, lm_stat, R = 5L, verbose = "yes"), "must be")
  expect_err(fwb(d, lm_stat, R = 5L, verbose = FALSE, simple = "yes"), "must be")
  expect_err(fwb(d, lm_stat, R = 5L, verbose = FALSE, wtype = "nope"),
             "should be one of")
  expect_err(fwb(d, lm_stat, R = 5L, verbose = FALSE, wtype = c("exp", "beta")),
             "must be")

  #`drop0` accepts `TRUE`, `FALSE`, and `NA`, but only for the weight types that can
  #produce a zero.
  expect_err(fwb(d, lm_stat, R = 5L, verbose = FALSE, wtype = "multinom",
                 drop0 = "yes"),
             "must be")

  for (dz in list(TRUE, FALSE, NA)) {
    expect_no_error({
      f <- fwb(d, lm_stat, R = 5L, verbose = FALSE, wtype = "poisson", drop0 = dz)
    })
  }

  #For the continuous types `drop0` is silently ignored rather than rejected, since
  #there is nothing for it to do.
  expect_no_error({
    f <- fwb(d, lm_stat, R = 5L, verbose = FALSE, wtype = "exp", drop0 = "nonsense")
  })
})

test_that("`drop0` controls what `statistic` receives", {
  d <- fwb_test_df()

  #A statistic that reports what it was handed, so the three `drop0` settings can be
  #distinguished without inferring anything from the estimates.
  probe <- function(data, w) {
    c(n = nrow(data), n_na = sum(is.na(w)), n_zero = sum(w == 0, na.rm = TRUE))
  }

  set.seed(95, "L'Ecuyer-CMRG")
  keep <- fwb(d, probe, R = 10L, verbose = FALSE, wtype = "multinom",
              drop0 = FALSE, simple = TRUE)

  expect_true(all(keep[["t"]][, "n"] == nrow(d)))
  expect_true(all(keep[["t"]][, "n_na"] == 0))
  expect_true(all(keep[["t"]][, "n_zero"] > 0))

  set.seed(95, "L'Ecuyer-CMRG")
  dropped <- fwb(d, probe, R = 10L, verbose = FALSE, wtype = "multinom",
                 drop0 = TRUE, simple = TRUE)

  expect_true(all(dropped[["t"]][, "n"] < nrow(d)))
  expect_true(all(dropped[["t"]][, "n_zero"] == 0))
  expect_true(all(dropped[["t"]][, "n_na"] == 0))

  set.seed(95, "L'Ecuyer-CMRG")
  na_out <- fwb(d, probe, R = 10L, verbose = FALSE, wtype = "multinom",
                drop0 = NA, simple = TRUE)

  expect_true(all(na_out[["t"]][, "n"] == nrow(d)))
  expect_true(all(na_out[["t"]][, "n_zero"] == 0))
  expect_equal(na_out[["t"]][, "n_na"], keep[["t"]][, "n_zero"],
               ignore_attr = TRUE)

  #All three routes must give the same estimates for a `statistic` that handles the
  #zeros correctly, which is the whole justification for offering the choice.
  set.seed(95, "L'Ecuyer-CMRG")
  e_keep <- fwb(d, lm_stat, R = 10L, verbose = FALSE, wtype = "multinom",
                drop0 = FALSE, simple = TRUE)

  set.seed(95, "L'Ecuyer-CMRG")
  e_drop <- fwb(d, lm_stat, R = 10L, verbose = FALSE, wtype = "multinom",
                drop0 = TRUE, simple = TRUE)

  set.seed(95, "L'Ecuyer-CMRG")
  e_na <- fwb(d, lm_stat, R = 10L, verbose = FALSE, wtype = "multinom",
              drop0 = NA, simple = TRUE)

  expect_same_t(e_drop, e_keep)
  expect_same_t(e_na, e_keep)
})

test_that("extra arguments reach `statistic`", {
  d <- fwb_test_df()

  extra_stat <- function(data, w, mult, add = 0) {
    mult * coef(lm(y ~ x1, data = data, weights = w)) + add
  }

  set.seed(96, "L'Ecuyer-CMRG")
  f <- fwb(d, extra_stat, R = 10L, verbose = FALSE, mult = 2, add = 1)

  set.seed(96, "L'Ecuyer-CMRG")
  ref <- fwb(d, function(data, w) coef(lm(y ~ x1, data = data, weights = w)),
             R = 10L, verbose = FALSE)

  expect_equal(unname(f[["t"]]), 2 * unname(ref[["t"]]) + 1,
               tolerance = fwb_eps())

  #A missing required argument surfaces through the `statistic` error wrapper.
  expect_err(fwb(d, extra_stat, R = 5L, verbose = FALSE),
             "there was a problem running the function supplied to")
})

test_that("unnamed statistics get default names", {
  d <- fwb_test_df()

  set.seed(97, "L'Ecuyer-CMRG")
  f <- fwb(d, function(data, w) unname(coef(lm(y ~ x1, data = data, weights = w))),
           R = 8L, verbose = FALSE)

  expect_identical(names(f[["t0"]]), c("t1", "t2"))
  expect_identical(colnames(f[["t"]]), c("t1", "t2"))
})

test_that("the returned `<fwb>` object has the documented structure", {
  d <- fwb_test_df()

  set.seed(98, "L'Ecuyer-CMRG")
  f <- fwb(d, lm_stat, R = 12L, verbose = FALSE)

  expect_identical(names(f),
                   c("t0", "t", "R", "data", "statistic", "call",
                     "cluster", "strata", "wtype"))
  expect_s3_class(f, c("fwb", "boot"))

  expect_identical(f[["data"]], d)
  expect_identical(f[["statistic"]], lm_stat)
  expect_identical(f[["R"]], 12L)
  expect_null(f[["cluster"]])
  expect_null(f[["strata"]])
  expect_identical(f[["wtype"]], "exp")

  #`call` is the call as written, which is what `print()` shows.
  expect_true(is.call(f[["call"]]))
  expect_identical(as.character(f[["call"]][[1L]]), "fwb")

  #The attributes that `fwb.array()` reads.
  expect_identical(attr(f, "boot_type", TRUE), "fwb")
  expect_true(attr(f, "simple", TRUE))

  #Whatever is needed to re-generate the weights goes in the `"seeds"` attribute. With
  #`simple = TRUE` that is one recorded L'Ecuyer stream per replicate, which is what makes
  #the weights recoverable without knowing anything about how the run was parallelized --
  #so the backend itself is deliberately *not* recorded.
  seeds <- attr(f, "seeds", TRUE)

  expect_true(is.matrix(seeds))
  expect_identical(nrow(seeds), 12L)
  expect_identical(typeof(seeds), "integer")
  expect_false(anyNA(seeds))
  expect_identical(anyDuplicated(seeds), 0L)
  expect_null(attr(f, "cl", TRUE))

  #`simple = FALSE` draws the whole matrix in one call, so a single generator state is
  #enough and `"seeds"` holds that one state rather than a matrix of them.
  set.seed(98, "L'Ecuyer-CMRG")
  fns <- fwb(d, lm_stat, R = 12L, verbose = FALSE, simple = FALSE)

  seed <- attr(fns, "seeds", TRUE)

  expect_type(seed, "integer")
  expect_false(is.matrix(seed))
})

test_that("`vcovFWB()` validates its arguments", {
  d <- fwb_test_df()
  fit <- lm(y ~ x1 + x2, data = d)

  expect_err(vcovFWB(fit, R = 2.5), "must be")
  expect_err(vcovFWB(fit, R = 10L, start = "yes"), "must be")
  expect_err(vcovFWB(fit, R = 10L, fix = "yes"), "must be")
  expect_err(vcovFWB(fit, R = 10L, verbose = "yes"), "must be")
  expect_err(vcovFWB(fit, R = 10L, use = 1), "must be a string")
  expect_err(vcovFWB(fit, R = 10L, .coef = "coef"), "must be a function")
  expect_err(vcovFWB(fit, R = 10L, wtype = "nope"), "should be one of")

  #A `cluster` of the wrong length.
  expect_err(vcovFWB(fit, R = 10L, cluster = seq_len(5L)),
             "do not match")

  #`NA`s in `cluster` cannot be handled here; the user has to resolve them.
  bad_cluster <- as.integer(d[["g"]])
  is.na(bad_cluster[1L]) <- TRUE

  expect_err(vcovFWB(fit, R = 10L, cluster = bad_cluster), "cannot handle")
})

test_that("`vcovFWB()` requires `.coef` to return a numeric vector", {
  d <- fwb_test_df()
  fit <- lm(y ~ x1 + x2, data = d)

  #The default message points at `.coef` as the way out, because the usual cause is
  #a model whose `coef()` returns a matrix (`nnet::multinom()`, and similar).
  expect_err(vcovFWB(fit, R = 10L, .coef = function(x) matrix(1:4, 2L)),
             "must return a numeric vector")

  #A custom `.coef` that does return a vector is accepted.
  expect_no_error({
    v <- vcovFWB(fit, R = 10L, .coef = function(x) stats::coef(x)[1:2])
  })

  expect_identical(dim(v), c(2L, 2L))
})

test_that("`vcovFWB()` warns when the model ignores the weights", {
  d <- fwb_test_df()

  #A model refit without the bootstrap weights produces identical coefficients every
  #replicate, so every variance is 0. That is silent nonsense unless it is flagged:
  #the returned matrix looks like an extremely precise estimate.
  unweighted <- lm(y ~ x1, data = d, weights = rep(0, nrow(d)) + 1)

  #`.coef` that discards the fit entirely is the cleanest way to force the
  #degenerate case without depending on a model that mishandles `weights`.
  set.seed(99, "L'Ecuyer-CMRG")
  expect_wrn(v <- vcovFWB(unweighted, R = 10L,
                          .coef = function(x) c(a = 1, b = 2)),
             "all variances and covariances are 0")

  expect_equal(v, matrix(0, 2L, 2L, dimnames = list(c("a", "b"), c("a", "b"))),
               tolerance = fwb_eps())
})

test_that("`vcovFWB()` accepts a `cluster` formula, vector, or data frame", {
  d <- fwb_test_df()
  fit <- lm(y ~ x1 + x2, data = d)

  set.seed(100, "L'Ecuyer-CMRG")
  by_formula <- vcovFWB(fit, cluster = ~ g, R = 20L)

  set.seed(100, "L'Ecuyer-CMRG")
  by_vector <- vcovFWB(fit, cluster = d[["g"]], R = 20L)

  set.seed(100, "L'Ecuyer-CMRG")
  by_df <- vcovFWB(fit, cluster = d["g"], R = 20L)

  expect_equal(by_formula, by_vector, tolerance = fwb_eps())
  expect_equal(by_formula, by_df, tolerance = fwb_eps())

  #Clustering must change the answer -- otherwise the argument is being dropped.
  set.seed(100, "L'Ecuyer-CMRG")
  unclustered <- vcovFWB(fit, R = 20L)

  expect_not_equal(by_formula, unclustered)

  #Multi-way clustering applies the inclusion-exclusion correction, so the result is
  #not simply either one-way matrix.
  set.seed(100, "L'Ecuyer-CMRG")
  two_way <- vcovFWB(fit, cluster = ~ g + s, R = 20L)

  expect_identical(dim(two_way), dim(by_formula))
  expect_not_equal(two_way, by_formula)
})

test_that("`vcovFWB()` matches `fwb()` on the same weights", {
  d <- fwb_test_df()

  #The two functions draw weights the same way, so given the same seed and the same
  #model they must produce the same covariance. This is the invariant that lets
  #`vcovFWB()` be described as a shortcut for `fwb()` plus `vcov()`.
  for (fit_spec in list(
    list(fit = lm(y ~ x1 + x2, data = d),
         stat = function(data, w) coef(lm(y ~ x1 + x2, data = data, weights = w))),
    list(fit = glm(yb ~ x1 + x2, data = d, family = quasibinomial()),
         stat = function(data, w) {
           coef(glm(yb ~ x1 + x2, data = data, family = quasibinomial(),
                    weights = w))
         }))) {

    for (wt in c("exp", "multinom", "poisson", "beta")) {
      set.seed(101, "L'Ecuyer-CMRG")
      v <- vcovFWB(fit_spec[["fit"]], R = 40L, wtype = wt)

      set.seed(101, "L'Ecuyer-CMRG")
      f <- fwb(d, fit_spec[["stat"]], R = 40L, verbose = FALSE, wtype = wt,
               simple = TRUE)

      expect_equal(v, vcov(f), tolerance = 1e-6,
                   info = paste(class(fit_spec[["fit"]])[1L], wt))
    }
  }
})

test_that("`vcovFWB()` multiplies existing model weights by the bootstrap weights", {
  d <- fwb_test_df()
  d[["sw"]] <- 1 + as.numeric(d[["x2"]])

  fit <- lm(y ~ x1, data = d, weights = sw)

  set.seed(102, "L'Ecuyer-CMRG")
  v <- vcovFWB(fit, R = 40L)

  set.seed(102, "L'Ecuyer-CMRG")
  f <- fwb(d, function(data, w) coef(lm(y ~ x1, data = data, weights = w * sw)),
           R = 40L, verbose = FALSE)

  expect_equal(v, vcov(f), tolerance = 1e-6)

  #Ignoring the sampling weights would give a different answer, so the comparison
  #above is informative.
  set.seed(102, "L'Ecuyer-CMRG")
  f_noSW <- fwb(d, function(data, w) coef(lm(y ~ x1, data = data, weights = w)),
                R = 40L, verbose = FALSE)

  expect_not_equal(v, vcov(f_noSW))
})

test_that("`vcovFWB()` requires `R > 1`", {
  #`R` used to be checked only for being a count, so `R = 0` ran zero replicates and
  #failed inside `stats::cov()` with "supply both 'x' and 'y' or a matrix-like 'x'".
  #The bound is 1 rather than 0 because a covariance needs at least two replicates --
  #`R = 1` would give a matrix of `NA`s.
  d <- fwb_test_df()
  fit <- lm(y ~ x1 + x2, data = d)

  expect_err(vcovFWB(fit, R = 0L), "must be greater than 1")
  expect_err(vcovFWB(fit, R = 1L), "must be greater than 1")

  set.seed(129, "L'Ecuyer-CMRG")
  expect_no_error({
    v <- vcovFWB(fit, R = 2L)
  })

  expect_true(all(is.finite(v)))
})

test_that("`fwb()` refuses a dataset with nothing to resample", {
  #One unit has exactly one possible bootstrap sample, so every replicate returns the
  #original estimate and every standard error is 0. Returning zeros silently is worse
  #than refusing.
  st <- function(data, w) c(m = w_mean(data[["x"]], w))

  expect_err(fwb(data.frame(x = 1), st, R = 9L, verbose = FALSE),
             "must have more than one unit")

  #Two is enough.
  expect_no_error(fwb(data.frame(x = c(1, 2)), st, R = 9L, wtype = "exp",
                      verbose = FALSE))

  #The same holds one level up for clusters: the cluster bootstrap resamples clusters.
  d <- data.frame(x = 1:6, g = factor(rep(1L, 6L)))

  expect_err(fwb(d, st, R = 9L, cluster = g, verbose = FALSE),
             "more than one cluster")

  expect_no_error({
    d[["g2"]] <- factor(rep(1:2, each = 3L))
    fwb(d, st, R = 9L, cluster = g2, wtype = "exp", verbose = FALSE)
  })
})

test_that("`vcovFWB()` refuses a model with nothing to resample", {
  fit <- lm(y ~ 1, data = data.frame(y = 3))

  expect_err(vcovFWB(fit, R = 9L), "more than one unit")

  #And a clustering variable with a single cluster, which is the same degeneracy.
  d <- data.frame(y = c(1, 3, 2, 5), x = c(1, 2, 3, 4), g = rep(1L, 4L))
  fit2 <- lm(y ~ x, data = d)

  expect_err(vcovFWB(fit2, cluster = ~g, R = 9L), "more than one cluster")
})

test_that("`statistic` is validated, but only when its formals allow it", {
  d <- fwb_test_df()

  #Cannot receive the weights at all.
  expect_err(fwb(d, function(data) 1, R = 5L, verbose = FALSE),
             "must accept at least two arguments")

  #An argument sharing a name with one of `fwb()`'s own can never be supplied: `fwb()`
  #matches it first.
  expect_err(fwb(d, function(data, w, verbose) 1, R = 5L, verbose = FALSE),
             "cannot have argument")

  expect_err(fwb(d, function(data, w, cl, wtype) 1, R = 5L, verbose = FALSE),
             "cannot have arguments")

  #The first two arguments are passed positionally, so their names are unconstrained --
  #even when they collide.
  expect_no_error(fwb(d, function(R, wtype) c(m = w_mean(R[["x1"]], wtype)),
                      R = 5L, verbose = FALSE))

  #Extra arguments that do not collide are fine and still reachable through `...`.
  expect_no_error({
    f <- fwb(d, function(data, w, k) c(m = k * w_mean(data[["x1"]], w)),
             R = 5L, verbose = FALSE, k = 2)
  })

  #A function with `...` is exempt from every check above, because it can accept
  #whatever it is handed and nothing about its formals can be shown to be wrong. This is
  #what lets other packages wrap `statistic` without `fwb` knowing they exist: any
  #wrapper has to forward arguments, so any wrapper has `...`.
  expect_no_error(fwb:::check_statistic(function(...) 1))
  expect_no_error(fwb:::check_statistic(function(..., verbose) 1))

  #Shaped like *progressify*'s wrapper, but the exemption is keyed on `...`, not on any
  #knowledge of that package.
  expect_no_error(fwb:::check_statistic(
    function(..., ...FUN, .progressr_progressor) 1
  ))
})
