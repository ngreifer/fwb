#`vcovFWB()`'s refitting machinery.
#
#`vcovFWB()` picks one of three ways to refit the model `R` times: `.lm.fit()` for a
#plain `lm`, the original `glm` fitting method for a `glm`, or a general `update()`
#call for anything else. The choice is made by refitting once and checking whether the
#model matrix changed -- a `w_*()` term in the formula makes it change, which forces
#the general path.
#
#`test-vcovFWB.R` checks that the covariances come out right. This file checks *which
#path was taken* and what happens at the edges of that decision, because the fast
#paths and the general path are separate implementations of the same thing and can
#drift apart.

#Which of the three refitting paths `make.bootfit()` chose, read off the body of the
#function it returns. `.env` has to be the frame where the model's `data` lives,
#because the probe refit re-evaluates the model call there.
bootfit_path <- function(fit, drop0 = FALSE, wtype = "exp",
                         .coef = stats::coef, .env = parent.frame()) {
  bootfit <- fwb:::make.bootfit(fit, factor(seq_len(fwb:::nobs0(fit))),
                                start = NULL, drop0 = drop0,
                                fwb:::make_gen_weights(wtype),
                                .coef, .env)

  body_text <- paste(deparse(body(bootfit)), collapse = " ")

  if (grepl(".lm.fit", body_text, fixed = TRUE)) "lm"
  else if (grepl("safe.glm.fit", body_text, fixed = TRUE)) "glm"
  else "update"
}

test_that("`lm` and `glm` models take the fast refitting paths", {
  #`make.bootfit()` draws a set of weights for its probe refit, so the stream has to
  #exist before it is called.
  set.seed(118, "L'Ecuyer-CMRG")

  d <- fwb_test_df()
  env <- environment()

  expect_identical(bootfit_path(lm(y ~ x1 + x2, data = d), .env = env), "lm")
  expect_identical(bootfit_path(glm(yb ~ x1 + x2, data = d,
                                    family = quasibinomial()), .env = env),
                   "glm")

  #A `w_*()` term changes the model matrix from replicate to replicate, so the fast
  #paths -- which precompute `x` and `y` once -- cannot be used.
  expect_identical(bootfit_path(lm(y ~ x1 + w_center(x2), data = d), .env = env),
                   "update")

  #Nor can they when `.coef` disagrees with `stats::coef()`, since the fast paths
  #return raw fitting-function coefficients rather than running `.coef`.
  expect_identical(bootfit_path(lm(y ~ x1 + x2, data = d),
                                .coef = function(x) stats::coef(x)[1:2],
                                .env = env),
                   "update")
})

test_that("the fast and general paths give the same covariance", {
  set.seed(118, "L'Ecuyer-CMRG")

  d <- fwb_test_df()
  env <- environment()

  #`special_coef` is decided by comparing `.coef(fit)` against `stats::coef(fit)`, so
  #a wrapper that returns the same values takes the fast path anyway -- the check is on
  #the results, not on function identity.
  same_coef <- function(x) stats::coef(x)

  expect_identical(bootfit_path(lm(y ~ x1 + x2, data = d), .coef = same_coef,
                                .env = env),
                   "lm")

  #Reversing the coefficients is a bijection, so it forces the general `update()` path
  #while still computing every coefficient. Undoing the permutation must then recover
  #the fast path's covariance exactly: any difference is a difference between the two
  #implementations, not between the estimators.
  rev_coef <- function(x) rev(stats::coef(x))

  for (fit in list(lm(y ~ x1 + x2, data = d),
                   glm(yb ~ x1 + x2, data = d, family = quasibinomial()))) {

    expect_identical(bootfit_path(fit, .coef = rev_coef, .env = env), "update")

    set.seed(118, "L'Ecuyer-CMRG")
    fast <- vcovFWB(fit, R = 40L)

    set.seed(118, "L'Ecuyer-CMRG")
    general <- vcovFWB(fit, R = 40L, .coef = rev_coef)

    k <- rev(seq_len(nrow(fast)))

    expect_equal(fast, general[k, k], tolerance = 1e-6, info = class(fit)[1L])
  }
})

test_that("`drop0` does not change the covariance for `lm` or `glm`", {
  d <- fwb_test_df()

  #Zeroing a unit out and dropping it are equivalent for a weighted least squares or IRLS
  #fit, so all three `drop0` settings must agree. This is now a test of the fast paths'
  #own `drop0` handling, which was unreachable until the probe stopped applying `drop0`
  #-- and which was wrong when it was written: `drop0 = NA` handed `NA` weights straight
  #to the fitting function, which cannot take them.
  for (fit in list(lm(y ~ x1 + x2, data = d),
                   glm(yb ~ x1 + x2, data = d, family = quasibinomial()))) {

    set.seed(119, "L'Ecuyer-CMRG")
    keep <- vcovFWB(fit, R = 40L, wtype = "multinom", drop0 = FALSE)

    set.seed(119, "L'Ecuyer-CMRG")
    dropped <- vcovFWB(fit, R = 40L, wtype = "multinom", drop0 = TRUE)

    set.seed(119, "L'Ecuyer-CMRG")
    na_out <- vcovFWB(fit, R = 40L, wtype = "multinom", drop0 = NA)

    expect_equal(dropped, keep, tolerance = 1e-6, info = class(fit)[1L])
    expect_equal(na_out, keep, tolerance = 1e-6, info = class(fit)[1L])
  }
})

test_that("`drop0` keeps the fast refitting paths", {
  #The probe that chooses between the fast paths and the general `update()` path used to
  #refit *with* `drop0` applied, and both `subset = w > 0` (`drop0 = TRUE`) and `NA`
  #weights (`drop0 = NA`) change the number of rows in the model matrix. The probe read
  #that as "the model matrix depends on the weights" and sent every model down the general
  #path, which cost the `.lm.fit()`/`glm.fit()` speedup on all `R` replicates and left the
  #fast paths' own `drop0` handling unreachable.
  set.seed(119, "L'Ecuyer-CMRG")

  d <- fwb_test_df()
  env <- environment()

  for (spec in list(list(fit = lm(y ~ x1 + x2, data = d), path = "lm"),
                    list(fit = glm(yb ~ x1 + x2, data = d,
                                   family = quasibinomial()), path = "glm"))) {
    for (dz in list(FALSE, TRUE, NA)) {
      expect_identical(bootfit_path(spec[["fit"]], drop0 = dz,
                                    wtype = "multinom", .env = env),
                       spec[["path"]],
                       info = paste(spec[["path"]], "drop0 =", dz))
    }
  }

  #A model whose model matrix genuinely does depend on the weights must still take the
  #general path, whatever `drop0` is -- otherwise the probe has been defanged rather
  #than fixed.
  for (dz in list(FALSE, TRUE, NA)) {
    expect_identical(bootfit_path(lm(y ~ x1 + w_center(x2), data = d), drop0 = dz,
                                  wtype = "multinom", .env = env),
                     "update",
                     info = paste("w_center, drop0 =", dz))
  }
})

test_that("the probe refit does not consume the RNG", {
  set.seed(120, "L'Ecuyer-CMRG")

  d <- fwb_test_df()
  fit <- lm(y ~ x1 + x2, data = d)

  #`make.bootfit()` draws one set of weights to compare model matrices. It wraps that
  #in `with_seed_preserved()`, so building the fitting function must leave the stream
  #exactly where it was -- otherwise `vcovFWB()` would not be reproducible from a
  #plain `set.seed()`.
  before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)

  invisible(bootfit_path(fit, .coef = function(x) stats::coef(x),
                         .env = environment()))

  expect_identical(get(".Random.seed", envir = globalenv(), inherits = FALSE),
                   before)
})

test_that("`start = TRUE` does not change the estimates", {
  d <- fwb_test_df()

  fit <- glm(yb ~ x1 + x2, data = d, family = quasibinomial())

  #`start` only supplies initial values to the IRLS iteration, so it may change the
  #speed but must not change the converged coefficients beyond tolerance.
  set.seed(121, "L'Ecuyer-CMRG")
  without <- vcovFWB(fit, R = 40L, start = FALSE)

  set.seed(121, "L'Ecuyer-CMRG")
  with_start <- vcovFWB(fit, R = 40L, start = TRUE)

  expect_equal(with_start, without, tolerance = 1e-5)

  #`start` is ignored for non-`glm` models rather than being passed on and erroring.
  lm_fit <- lm(y ~ x1 + x2, data = d)

  set.seed(121, "L'Ecuyer-CMRG")
  expect_no_error({
    v <- vcovFWB(lm_fit, R = 20L, start = TRUE)
  })
})

test_that("`vcovFWB()` refits models that need `update()`", {
  skip_if_not_installed("survival")

  d <- fwb_test_df()
  d[["time"]] <- exp(d[["y"]])

  fit <- survival::coxph(survival::Surv(time, yb) ~ x1 + x2, data = d)

  set.seed(122, "L'Ecuyer-CMRG")
  expect_no_error({
    v <- vcovFWB(fit, R = 30L)
  })

  expect_identical(dimnames(v), list(names(coef(fit)), names(coef(fit))))

  set.seed(122, "L'Ecuyer-CMRG")
  f <- fwb(d, function(data, w) {
    coef(survival::coxph(survival::Surv(time, yb) ~ x1 + x2, data = data,
                         weights = w))
  }, R = 30L, verbose = FALSE)

  expect_equal(v, vcov(f), tolerance = 1e-6)
})

test_that("`.coef` handles models whose `coef()` is not a vector", {
  skip_if_not_installed("nnet")

  d <- fwb_test_df()
  d[["grp"]] <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))

  fit <- nnet::multinom(grp ~ x1, data = d, trace = FALSE)

  #`coef()` on a multinomial fit returns a matrix, which `vcovFWB()` cannot use --
  #and the error is supposed to point at `.coef` as the fix.
  expect_err(vcovFWB(fit, R = 10L), "see the")

  coef_multinom <- function(x) {
    p <- t(coef(x))

    setNames(as.vector(p),
             paste(colnames(p)[col(p)], rownames(p)[row(p)], sep = ":"))
  }

  set.seed(123, "L'Ecuyer-CMRG")
  expect_no_error({
    v <- vcovFWB(fit, R = 20L, .coef = coef_multinom)
  })

  expect_identical(dim(v), c(length(coef_multinom(fit)),
                             length(coef_multinom(fit))))
  expect_identical(colnames(v), names(coef_multinom(fit)))
})

test_that("`fix = TRUE` projects a non-PSD covariance to the nearest PSD one", {
  #Multi-way clustering combines one-way covariances with alternating signs, so the
  #result need not be positive semi-definite -- especially at small `R`. `fix = TRUE`
  #exists for exactly that case.
  #
  #It used to keep only `eigen(...)$values` and then reach for `eig$values`/`eig$vectors`
  #on that numeric vector, so any call that actually reached the repair failed with
  #"$ operator is invalid for atomic vectors" -- a no-op when the matrix was already PSD
  #and an error exactly when it was supposed to act.
  n <- 300L

  dd <- data.frame(x = qnorm(ppoints(n)),
                   y = sinpi(seq_len(n) / 11),
                   f1 = factor(rep(seq_len(6L), length.out = n)),
                   f2 = factor(rep(seq_len(7L), length.out = n)))

  m <- lm(y ~ x, data = dd)

  set.seed(124, "L'Ecuyer-CMRG")
  raw <- vcovFWB(m, cluster = ~ f1 + f2, R = 10L, fix = FALSE)

  #The setup has to actually produce a non-PSD matrix, or the test proves nothing.
  expect_true(any(eigen(raw, symmetric = TRUE)$values < 0))

  set.seed(124, "L'Ecuyer-CMRG")
  expect_no_error({
    fixed <- vcovFWB(m, cluster = ~ f1 + f2, R = 10L, fix = TRUE)
  })

  expect_true(all(eigen(fixed, symmetric = TRUE)$values >= -fwb_eps()))
  expect_identical(dimnames(fixed), dimnames(raw))
  expect_true(isSymmetric(fixed))
})

test_that("`fix = TRUE` is a no-op on a positive semi-definite covariance", {
  d <- fwb_test_df()
  fit <- lm(y ~ x1 + x2, data = d)

  set.seed(125, "L'Ecuyer-CMRG")
  unfixed <- vcovFWB(fit, R = 40L, fix = FALSE)

  expect_true(all(eigen(unfixed, symmetric = TRUE)$values >= 0))

  set.seed(125, "L'Ecuyer-CMRG")
  fixed <- vcovFWB(fit, R = 40L, fix = TRUE)

  expect_equal(fixed, unfixed, tolerance = fwb_eps())
})

test_that("`vcovFWB()` handles models fit on incomplete data", {
  #Weights are computed at the length of the rows the fit used, but `update()`
  #re-evaluates the model call against the full data frame. When `na.action` dropped rows
  #the two lengths differed and the refit failed with "variable lengths differ (found for
  #'(weights)')", so `vcovFWB()` could not be used at all on a model fit on data with
  #missing values -- with or without `cluster`, for `lm` and `glm` alike.
  #
  #`fwb()` was never affected: `statistic` decides for itself what to do with `NA`s.
  d <- fwb_test_df()
  is.na(d[["x1"]][c(2L, 20L, 45L)]) <- TRUE

  fits <- list(lm = lm(y ~ x1 + x2, data = d),
               glm = glm(yb ~ x1 + x2, data = d, family = quasibinomial()),
               `lm + w_center` = lm(y ~ x1 + w_center(x2), data = d))

  for (nm in names(fits)) {
    fit <- fits[[nm]]

    expect_identical(fwb:::nobs0(fit), nrow(d) - 3L, info = nm)

    set.seed(126, "L'Ecuyer-CMRG")
    expect_no_error({
      v0 <- vcovFWB(fit, R = 30L)
    })

    expect_identical(dim(v0), c(3L, 3L), info = nm)
    expect_true(all(is.finite(v0)), info = nm)

    #`cluster` is supplied at full length; `vcovFWB()` drops the rows the model dropped,
    #using the fit's own `na.action` record -- the same record used to put the weights
    #back where they came from.
    set.seed(126, "L'Ecuyer-CMRG")
    expect_no_error({
      v <- vcovFWB(fit, cluster = ~ g, R = 30L)
    })

    expect_identical(dim(v), c(3L, 3L), info = nm)
  }

  #Missingness in the response, not just a predictor.
  d2 <- fwb_test_df()
  is.na(d2[["y"]][c(5L, 30L)]) <- TRUE

  set.seed(126, "L'Ecuyer-CMRG")
  expect_no_error({
    v <- vcovFWB(lm(y ~ x1 + x2, data = d2), R = 30L)
  })

  expect_true(all(is.finite(v)))

  #And the answer is the one you would get by dropping the incomplete rows yourself,
  #which is what makes the padding correct rather than merely non-erroring.
  complete <- d[stats::complete.cases(d[c("y", "x1", "x2")]), ]

  set.seed(101, "L'Ecuyer-CMRG")
  v <- vcovFWB(fits[["lm"]], R = 40L)

  set.seed(101, "L'Ecuyer-CMRG")
  f <- fwb(complete, lm_stat, R = 40L, verbose = FALSE)

  expect_equal(v, vcov(f), tolerance = 1e-6)
})

test_that("`vcovFWB()` muffles the per-replicate non-integer-weight warnings", {
  d <- fwb_test_df()

  #Fractional weights in a binomial `glm()` produce "non-integer #successes in a
  #binomial glm!" on every replicate. Relaying `R` copies of that would bury anything
  #real, so `safe.glm.fit()` swallows it -- but only it, and only for this reason.
  fit <- glm(yb ~ x1 + x2, data = d, family = binomial())

  warnings_seen <- character()

  set.seed(127, "L'Ecuyer-CMRG")
  v <- withCallingHandlers(vcovFWB(fit, R = 20L),
                           warning = function(w) {
                             warnings_seen <<- c(warnings_seen,
                                                 conditionMessage(w))
                             invokeRestart("muffleWarning")
                           })

  expect_identical(dim(v), c(3L, 3L))

  #None gets through. The probe refit in `make.bootfit()` calls `update()` directly
  #rather than `safe.glm.fit()`, so it needs its own suppression -- without it exactly
  #one warning leaked from every `vcovFWB()` call, which reads like a real diagnostic.
  expect_identical(sum(grepl("non-integer", warnings_seen, fixed = TRUE)), 0L)
})

test_that("the fast `lm` path handles a rank-deficient replicate like `lm()` does", {
  #`.lm.fit()` is the raw decomposition. Unlike `lm.fit()` it leaves the coefficients in
  #pivoted order and leaves the aliased ones at whatever the decomposition produced
  #(usually 0), so a replicate whose weights leave the design rank-deficient needs both
  #undone. `x2` is 1 for a single unit, so any resample that drops that unit makes the
  #column all-zero and therefore aliased; putting it *before* `x1` in the formula means
  #the pivot reorders, which is the damaging case.
  d <- data.frame(x1 = c(1, 2, 3, 4, 5, 6),
                  x2 = c(1, 0, 0, 0, 0, 0),
                  y  = c(3, 1, 4, 2, 6, 5))

  fit <- lm(y ~ x2 + x1, data = d)

  x <- model.matrix(fit)
  yv <- d[["y"]]

  fast <- function(w) {
    s <- sqrt(w)
    z <- .lm.fit(x * s, yv * s)
    cf <- z[["coefficients"]]

    if (z[["rank"]] < ncol(x)) {
      cf[seq.int(z[["rank"]] + 1L, ncol(x))] <- NA_real_
    }

    if (z[["pivoted"]]) {
      cf[z[["pivot"]]] <- cf
    }

    cf
  }

  set.seed(300, "L'Ecuyer-CMRG")
  gen <- fwb:::make_gen_weights("multinom")

  W <- t(vapply(seq_len(60L), function(i) drop(gen(6L, 1L)), numeric(6L)))

  #Only worth asserting if the situation actually arises in the sample of weights.
  deficient <- apply(W, 1L, function(w) .lm.fit(x * sqrt(w), yv * sqrt(w))[["rank"]] < 3L)
  reordered <- apply(W, 1L, function(w) .lm.fit(x * sqrt(w), yv * sqrt(w))[["pivoted"]])

  expect_true(any(deficient))
  expect_true(any(reordered))

  for (i in which(deficient)) {
    expect_equal(unname(fast(W[i, ])),
                 unname(coef(lm(y ~ x2 + x1, data = d, weights = W[i, ]))),
                 tolerance = fwb_eps(),
                 info = sprintf("replicate %d", i))
  }

  #The un-reordered raw result disagrees, which is what makes this worth a test: it
  #reports `x1`'s estimate against `x2`.
  i <- which(reordered)[1L]
  raw <- .lm.fit(x * sqrt(W[i, ]), yv * sqrt(W[i, ]))[["coefficients"]]

  expect_false(isTRUE(all.equal(unname(raw), unname(fast(W[i, ])))))
})

test_that("`vcovFWB()` matches `fwb()` when replicates are rank-deficient", {
  #The end-to-end consequence: an aliased coefficient left at 0 rather than `NA` cannot
  #be dropped by `use = "pairwise.complete.obs"`, so a replicate that estimated nothing
  #would be averaged in as though it had. `fwb()` and the general `update()` path both
  #give `NA`, so agreement with `fwb()` is the check.
  d <- data.frame(x1 = c(1, 2, 3, 4, 5, 6),
                  x2 = c(1, 0, 0, 0, 0, 0),
                  y  = c(3, 1, 4, 2, 6, 5))

  fit <- lm(y ~ x2 + x1, data = d)

  set.seed(301, "L'Ecuyer-CMRG")
  v <- vcovFWB(fit, R = 100L, wtype = "multinom")

  set.seed(301, "L'Ecuyer-CMRG")
  f <- suppressWarnings(
    fwb(d, function(data, w) coef(lm(y ~ x2 + x1, data = data, weights = w)),
        R = 100L, wtype = "multinom", simple = TRUE, verbose = FALSE)
  )

  expect_true(anyNA(f[["t"]]))

  expect_equal(unname(v),
               unname(cov(f[["t"]], use = "pairwise.complete.obs")),
               tolerance = fwb_eps())
})
