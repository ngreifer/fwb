#Edge cases for the `w_*()` functions.
#
#`test-w_mean.R` covers the main contract: agreement with the base R equivalents, and
#automatic weight capture inside `fwb()`. This file adds the cases where that contract
#has sharp edges -- missing values, zero weights, weights that reach the functions
#from `vcovFWB()` rather than `fwb()`, and the length changes that follow.

test_that("`w_*()` return unweighted values outside `fwb()`", {
  d <- fwb_test_df()
  x <- d[["y"]]

  #No weights anywhere means the base R answer, exactly.
  expect_equal(w_mean(x), mean(x), tolerance = fwb_eps())
  expect_equal(w_var(x), var(x), tolerance = fwb_eps())
  expect_equal(w_sd(x), sd(x), tolerance = fwb_eps())
  expect_equal(w_median(x), median(x), tolerance = fwb_eps())
  expect_equal(w_quantile(x), quantile(x), tolerance = fwb_eps())
  expect_equal(w_cov(cbind(x, d[["x1"]])), cov(cbind(x, d[["x1"]])),
               tolerance = fwb_eps())
  expect_equal(w_cor(cbind(x, d[["x1"]])), cor(cbind(x, d[["x1"]])),
               tolerance = fwb_eps())

  #Uniform weights give the same answer as no weights.
  w1 <- rep(1, length(x))

  expect_equal(w_mean(x, w1), mean(x), tolerance = fwb_eps())
  expect_equal(w_var(x, w1), var(x), tolerance = fwb_eps())
  expect_equal(w_median(x, w1), median(x), tolerance = fwb_eps())

  #And so do weights that are uniform up to a scale factor: the weights are
  #normalized to sum to 1 first.
  expect_equal(w_mean(x, w1 * 7), mean(x), tolerance = fwb_eps())
  expect_equal(w_var(x, w1 * 7), var(x), tolerance = fwb_eps())
})

test_that("`w_*()` are invariant to rescaling the weights", {
  d <- fwb_test_df()
  x <- d[["y"]]

  set.seed(107, "L'Ecuyer-CMRG")
  w <- rexp(length(x))

  #Every `w_*()` normalizes the weights to sum to 1, so the answer depends on the
  #weights only up to a positive scale factor.
  for (k in c(0.01, 3, 1000)) {
    expect_equal(w_mean(x, w * k), w_mean(x, w), tolerance = fwb_eps())
    expect_equal(w_var(x, w * k), w_var(x, w), tolerance = fwb_eps())
    expect_equal(w_median(x, w * k), w_median(x, w), tolerance = fwb_eps())
    expect_equal(w_quantile(x, w * k, probs = c(.1, .9)),
                 w_quantile(x, w, probs = c(.1, .9)), tolerance = fwb_eps())
  }
})

test_that("`na.rm` drops incomplete pairs", {
  d <- fwb_test_df()

  x <- d[["y"]]
  is.na(x[c(3L, 17L)]) <- TRUE

  set.seed(108, "L'Ecuyer-CMRG")
  w <- rexp(length(x))
  is.na(w[c(5L, 17L)]) <- TRUE

  #Without `na.rm` the missing values propagate, as they do in `mean()`/`var()`.
  expect_true(is.na(w_mean(x, w)))
  expect_true(is.na(w_var(x, w)))
  expect_true(is.na(w_sd(x, w)))

  #With `na.rm`, only rows complete in *both* `x` and `w` are used.
  keep <- !is.na(x) & !is.na(w)

  expect_equal(w_mean(x, w, na.rm = TRUE), w_mean(x[keep], w[keep]),
               tolerance = fwb_eps())
  expect_equal(w_var(x, w, na.rm = TRUE), w_var(x[keep], w[keep]),
               tolerance = fwb_eps())
  expect_equal(w_median(x, w, na.rm = TRUE), w_median(x[keep], w[keep]),
               tolerance = fwb_eps())
  expect_equal(w_quantile(x, w, probs = .25, na.rm = TRUE, names = FALSE),
               w_quantile(x[keep], w[keep], probs = .25, names = FALSE),
               tolerance = fwb_eps())
})

test_that("`w_cov()` handles missing values pairwise", {
  d <- fwb_test_df()

  X <- as.matrix(d[, c("y", "x1", "x2")])
  is.na(X[3L, 1L]) <- TRUE
  is.na(X[9L, 2L]) <- TRUE

  set.seed(109, "L'Ecuyer-CMRG")
  w <- rexp(nrow(X))

  #Without `na.rm`, the affected entries are missing.
  m <- w_cov(X, w)

  expect_true(anyNA(m))

  #With `na.rm`, each entry uses the rows complete for that pair -- so the (x2, x2)
  #variance, whose column has no missing values, is unaffected.
  m <- w_cov(X, w, na.rm = TRUE)

  expect_false(anyNA(m))
  expect_true(isSymmetric(m))
  expect_identical(dimnames(m), list(colnames(X), colnames(X)))

  expect_equal(m["x2", "x2"], w_var(X[, "x2"], w), tolerance = fwb_eps())

  keep <- !is.na(X[, "y"])

  expect_equal(m["y", "y"], w_var(X[keep, "y"], w[keep]), tolerance = fwb_eps())
})

test_that("`w_quantile()` interpolates the weighted CDF", {
  #With fewer than two positive weights there is nothing to interpolate between, so
  #`w_quantile()` documents a fallback to the unweighted `quantile()` rather than
  #returning the one positive-weight value for every `probs`.
  x <- c(1, 2, 3)

  expect_equal(w_quantile(x, c(0, 1, 0), probs = c(0, .5, 1), names = FALSE),
               quantile(x, probs = c(0, .5, 1), names = FALSE),
               tolerance = fwb_eps())

  #Weights of 0 drop their observations rather than distorting the interpolation, so
  #the result matches passing only the positive-weight rows.
  set.seed(110, "L'Ecuyer-CMRG")
  y <- rnorm(40L)
  w <- c(rep(0, 10L), rexp(30L))

  expect_equal(w_quantile(y, w, probs = c(.1, .5, .9)),
               w_quantile(y[w > 0], w[w > 0], probs = c(.1, .5, .9)),
               tolerance = fwb_eps())

  #Extremes are clamped to the range of `x` rather than extrapolated.
  q <- w_quantile(y, w, probs = c(0, 1), names = FALSE)

  expect_equal(q, range(y[w > 0]), tolerance = fwb_eps())

  #Unsorted input is handled.
  ord <- order(y)

  expect_equal(w_quantile(y, w, probs = .3), w_quantile(y[ord], w[ord], probs = .3),
               tolerance = fwb_eps())

  #Names follow `quantile()`'s convention unless suppressed.
  expect_named(w_quantile(y, w, probs = c(.025, .5)), c("2.5%", "50%"))
  expect_null(names(w_quantile(y, w, probs = .5, names = FALSE)))

  #With one positive weight there is nothing to interpolate, so it falls back to the
  #unweighted `quantile()`.
  one <- c(1, rep(0, 39L))

  expect_equal(w_quantile(y, one), quantile(y), tolerance = fwb_eps())

  expect_err(w_quantile(y, w, probs = 1.5), "must be between 0 and 1")
  expect_err(w_quantile(y, w, type = 6L), "must be equal to")
})

test_that("`w_median()` agrees with `w_quantile(probs = .5)`", {
  set.seed(111, "L'Ecuyer-CMRG")
  y <- rnorm(50L)
  w <- rexp(50L)

  expect_equal(w_median(y, w), w_quantile(y, w, probs = .5, names = FALSE),
               tolerance = fwb_eps())
})

test_that("`w_std()`, `w_scale()`, and `w_center()` compose consistently", {
  set.seed(112, "L'Ecuyer-CMRG")
  y <- rnorm(50L)
  w <- rexp(50L)

  expect_equal(w_std(y, w), (y - w_mean(y, w)) / w_sd(y, w),
               tolerance = fwb_eps())
  expect_equal(w_center(y, w), y - w_mean(y, w), tolerance = fwb_eps())

  #The standardized variable has weighted mean 0 and weighted variance 1.
  z <- w_std(y, w)

  expect_equal(w_mean(z, w), 0, tolerance = fwb_eps())
  expect_equal(w_var(z, w), 1, tolerance = fwb_eps())

  #`w_scale()` divides by the weighted standard deviation -- the *centered* second
  #moment -- matching both `?w_mean` and the unweighted branch. It used to use the
  #uncentered moment when weights were supplied.
  expect_equal(w_scale(y, w), y / w_sd(y, w), tolerance = fwb_eps())
})

test_that("`w_scale()` treats unit weights like no weights", {
  #`w_std(center = FALSE)` used to compute its scale factor two different ways: with no
  #weights `var(x)`, which centers, and with weights the uncentered
  #`sum(w * x^2) / (1 - sum(w^2))`. So `w_scale(x)` and `w_scale(x, rep(1, n))`
  #disagreed, unlike every other `w_*()` function, and inside `fwb()` -- where the
  #weights are never unit -- the formula in `?w_mean` was never the one applied.
  #
  #Note this deliberately does *not* match `scale(x, center = FALSE)`, which divides by
  #the root mean square. `w_scale()` follows its own documentation instead.
  set.seed(117, "L'Ecuyer-CMRG")
  y <- rnorm(50L)
  w1 <- rep(1, length(y))

  expect_equal(w_scale(y, w1), w_scale(y), tolerance = fwb_eps())
  expect_equal(w_scale(y, w1), y / sd(y), tolerance = fwb_eps())

  #The other two transformations already satisfy this.
  expect_equal(w_std(y, w1), w_std(y), tolerance = fwb_eps())
  expect_equal(w_center(y, w1), w_center(y), tolerance = fwb_eps())
})

test_that("`w_std()` keeps its input length when weights are 0 or missing", {
  #`w_std()` drops zero-weight and (with `na.rm = TRUE`, its default) missing rows when
  #computing the mean and scale, but must still return one value per element of `x`.
  #Anything else is fatal inside a model formula: `?w_mean` advertises
  #`lm(y ~ treat * w_center(x))` with `vcovFWB()`, which used to error with "variable
  #lengths differ" as soon as a weight was 0 -- that is, for `wtype = "multinom"` or
  #`"poisson"`.
  set.seed(113, "L'Ecuyer-CMRG")

  y <- rnorm(50L)
  w <- c(rep(0, 10L), rexp(40L))

  for (f in list(w_std, w_scale, w_center)) {
    expect_length(f(y, w), length(y))
  }

  #Missing values in `x` come back as missing rather than being dropped.
  yn <- y
  is.na(yn[3L]) <- TRUE

  z <- w_std(yn, rexp(50L))

  expect_length(z, length(yn))
  expect_true(is.na(z[3L]))
  expect_false(anyNA(z[-3L]))

  #The scale factor ignores the zero-weight rows, so the transformation is the one
  #computed from the positive-weight rows but applied to all of `x`.
  keep <- w > 0

  expect_equal(w_center(y, w), y - w_mean(y[keep], w[keep]), tolerance = fwb_eps())

  #And the whole point: a `w_*()` term in a formula now survives zero weights.
  d <- fwb_test_df()

  fit <- lm(y ~ x1 + w_center(x2), data = d)

  set.seed(113, "L'Ecuyer-CMRG")
  expect_no_error({
    v <- vcovFWB(fit, R = 20L, wtype = "multinom")
  })

  expect_identical(dim(v), c(3L, 3L))
})

test_that("`w_*()` capture the weights supplied by `vcovFWB()`", {
  d <- fwb_test_df()

  #The `w_*()` functions find the bootstrap weights through an option holding the
  #environment of the frame that called `statistic` (or that refit the model). This
  #test is the `vcovFWB()` half of that mechanism; `test-w_mean.R` covers the
  #`fwb()` half.
  #
  #A `w_center()` term forces `vcovFWB()` off its fast `.lm.fit()` path and onto the
  #general `update()` path, because the model matrix changes with the weights.
  fit <- lm(y ~ x1 * w_center(x2), data = d)

  set.seed(114, "L'Ecuyer-CMRG")
  v <- vcovFWB(fit, R = 40L)

  set.seed(114, "L'Ecuyer-CMRG")
  f <- fwb(d, function(data, w) {
    coef(lm(y ~ x1 * w_center(x2), data = data, weights = w))
  }, R = 40L, verbose = FALSE)

  expect_equal(v, vcov(f), tolerance = 1e-6)

  #Ignoring the weights inside `w_center()` gives a different answer, so the
  #comparison above is informative.
  set.seed(114, "L'Ecuyer-CMRG")
  f_fixed <- fwb(d, function(data, w) {
    coef(lm(y ~ x1 * I(x2 - mean(x2)), data = data, weights = w))
  }, R = 40L, verbose = FALSE)

  expect_not_equal(v, vcov(f_fixed))
})

test_that("`w_*()` do not leak weights outside the bootstrap", {
  d <- fwb_test_df()

  #The option that carries the weights is set with `rlang::local_options()` inside
  #the frame that calls `statistic`, so it must be gone once `fwb()` returns. If it
  #leaked, `w_mean()` at the console would silently report a weighted mean from the
  #last bootstrap replicate.
  set.seed(115, "L'Ecuyer-CMRG")
  invisible(fwb(d, function(data, w) c(m = w_mean(data[["y"]])), R = 8L,
                verbose = FALSE))

  expect_null(getOption("fwb_internal_w_env"))
  expect_equal(w_mean(d[["y"]]), mean(d[["y"]]), tolerance = fwb_eps())

  #Same after an error inside `statistic`.
  expect_err(fwb(d, function(data, w) stop("boom"), R = 8L, verbose = FALSE))

  expect_null(getOption("fwb_internal_w_env"))

  #And after `vcovFWB()`.
  set.seed(115, "L'Ecuyer-CMRG")
  invisible(vcovFWB(lm(y ~ x1, data = d), R = 8L))

  expect_null(getOption("fwb_internal_w_env"))
})

test_that("an explicit `w` overrides the captured bootstrap weights", {
  d <- fwb_test_df()

  w1 <- rep(1, nrow(d))

  set.seed(116, "L'Ecuyer-CMRG")
  f <- fwb(d, function(data, w) {
    c(auto = w_mean(data[["y"]]),
      unit = w_mean(data[["y"]], w1),
      explicit = w_mean(data[["y"]], w))
  }, R = 20L, verbose = FALSE)

  #The unit-weighted mean is constant across replicates; the other two vary and
  #agree with each other.
  expect_equal(sd(f[["t"]][, "unit"]), 0, tolerance = fwb_eps())
  expect_gt(sd(f[["t"]][, "auto"]), 0)
  expect_equal(f[["t"]][, "auto"], f[["t"]][, "explicit"], tolerance = fwb_eps())
  expect_equal(unname(f[["t"]][1L, "unit"]), mean(d[["y"]]), tolerance = fwb_eps())
})

test_that("`w_*()` reject unusable input", {
  expect_err(w_mean(), "must be supplied")
  expect_err(w_var(), "must be supplied")
  expect_err(w_sd(), "must be supplied")
  expect_err(w_median(), "must be supplied")
  expect_err(w_quantile(), "must be supplied")
  expect_err(w_cov(), "must be supplied")
  expect_err(w_cor(), "must be supplied")
  expect_err(w_std(), "must be supplied")

  x <- rnorm(10L)

  expect_err(w_mean(x, na.rm = "yes"), "must be")
  expect_err(w_quantile(list(1, 2)), "must be")
  expect_err(w_std(x, scale = "yes"), "must be")
  expect_err(w_std(x, center = "yes"), "must be")
})
