#Stratified and clustered bootstrapping.
#
#`?fwb` describes stratification as "essentially a separate bootstrap within each
#level of `strata`", which has a concrete numeric consequence: the weights sum to the
#stratum size within every stratum, in every replicate. The tests below check that
#through `fwb()` rather than through the generators directly (`test-wtypes.R` covers
#the generators), because the wiring between `fwb()` and the generator -- factor
#conversion, cluster/stratum collapsing -- is its own source of mistakes.

stratum_sums <- function(w, strata) {
  vapply(levels(strata), function(s) {
    rowSums(w[, strata == s, drop = FALSE])
  }, numeric(nrow(w)))
}

test_that("stratified weights sum to the stratum size in every replicate", {
  d <- fwb_test_df()
  sizes <- table(d[["s"]])

  for (wt in c("exp", "multinom", "mammen", "beta", "power")) {
    set.seed(72, "L'Ecuyer-CMRG")
    f <- fwb(d, w_identity_stat, R = 20L, verbose = FALSE, wtype = wt,
             simple = FALSE, strata = s)

    sums <- stratum_sums(f[["t"]], d[["s"]])

    expect_equal(sums,
                 matrix(as.numeric(sizes), nrow = 20L, ncol = length(sizes),
                        byrow = TRUE, dimnames = list(NULL, names(sizes))),
                 tolerance = fwb_eps(), info = wt)
  }
})

test_that("`strata` is stored and can be found in `data` or the caller", {
  d <- fwb_test_df()

  set.seed(73, "L'Ecuyer-CMRG")
  f <- fwb(d, lm_stat, R = 8L, verbose = FALSE, simple = FALSE, strata = s)

  #`strata` is stored as passed, not as the internal factor.
  expect_identical(f[["strata"]], d[["s"]])

  #`strata` is evaluated in `data` first, then in the calling frame.
  outside <- d[["s"]]

  set.seed(73, "L'Ecuyer-CMRG")
  f2 <- fwb(d, lm_stat, R = 8L, verbose = FALSE, simple = FALSE, strata = outside)

  expect_same_t(f2, f)
})

test_that("stratified bootstrapping works with `simple = TRUE`", {
  #`simple = TRUE` is the default for every `wtype` except `"multinom"`, so
  #`fwb(data, statistic, strata = s)` takes this path by default. It used to error,
  #because the generators indexed `w[, in_s]` without `drop = FALSE` and with `R = 1`
  #that drops the matrix to a vector before `rowMeans()`.
  d <- fwb_test_df()
  sizes <- table(d[["s"]])

  for (wt in c("exp", "mammen", "beta", "power")) {
    set.seed(74, "L'Ecuyer-CMRG")
    expect_no_error({
      f <- fwb(d, w_identity_stat, R = 12L, verbose = FALSE, wtype = wt,
               simple = TRUE, strata = s)
    })

    sums <- stratum_sums(f[["t"]], d[["s"]])

    expect_equal(as.vector(sums[1L, ]), as.numeric(sizes),
                 tolerance = fwb_eps(), info = wt)
  }
})

test_that("stratified bootstrapping works with `simple = TRUE` for the unaffected types", {
  d <- fwb_test_df()
  sizes <- table(d[["s"]])

  for (wt in c("poisson", "multinom")) {
    set.seed(74, "L'Ecuyer-CMRG")
    expect_no_error({
      f <- fwb(d, w_identity_stat, R = 12L, verbose = FALSE, wtype = wt,
               simple = TRUE, strata = s)
    })
  }

  #`"multinom"` still respects the stratum sizes when drawn one replicate at a time.
  set.seed(74, "L'Ecuyer-CMRG")
  f <- fwb(d, w_identity_stat, R = 12L, verbose = FALSE, wtype = "multinom",
           simple = TRUE, strata = s)

  sums <- stratum_sums(f[["t"]], d[["s"]])

  expect_equal(as.vector(sums[1L, ]), as.numeric(sizes), tolerance = fwb_eps())
})

test_that("a single-level `strata` is a no-op", {
  d <- fwb_test_df()
  d[["one"]] <- factor(rep("all", nrow(d)))

  set.seed(75, "L'Ecuyer-CMRG")
  f <- fwb(d, lm_stat, R = 12L, verbose = FALSE, simple = FALSE, strata = one)

  set.seed(75, "L'Ecuyer-CMRG")
  ref <- fwb(d, lm_stat, R = 12L, verbose = FALSE, simple = FALSE)

  expect_same_t(f, ref)
})

test_that("`\"poisson\"` ignores `strata` in `fwb()`", {
  d <- fwb_test_df()

  set.seed(76, "L'Ecuyer-CMRG")
  with_s <- fwb(d, w_identity_stat, R = 12L, verbose = FALSE, wtype = "poisson",
                simple = FALSE, strata = s)

  set.seed(76, "L'Ecuyer-CMRG")
  without <- fwb(d, w_identity_stat, R = 12L, verbose = FALSE, wtype = "poisson",
                 simple = FALSE)

  expect_same_t(with_s, without)
})

test_that("cluster weights are constant within cluster", {
  d <- fwb_test_df()

  for (simple in c(TRUE, FALSE)) {
    set.seed(77, "L'Ecuyer-CMRG")
    f <- fwb(d, w_identity_stat, R = 12L, verbose = FALSE, cluster = g,
             simple = simple)

    expect_true(all(apply(f[["t"]], 1L, function(w) {
      all(tapply(w, d[["g"]], function(z) fwb:::all_the_same(z)))
    })), info = paste("simple =", simple))

    #One distinct weight per cluster, not per unit.
    expect_true(all(apply(f[["t"]], 1L, function(w) {
      length(unique(round(w, 10L))) <= nlevels(d[["g"]])
    })))
  }

  expect_identical(nlevels(factor(d[["g"]])), 12L)
})

test_that("`cluster` is stored as a factor and can come from `data` or the caller", {
  d <- fwb_test_df()

  set.seed(78, "L'Ecuyer-CMRG")
  f <- fwb(d, lm_stat, R = 8L, verbose = FALSE, cluster = g)

  expect_s3_class(f[["cluster"]], "factor")
  expect_length(f[["cluster"]], nrow(d))
  expect_identical(levels(f[["cluster"]]), levels(d[["g"]]))

  outside <- d[["g"]]

  set.seed(78, "L'Ecuyer-CMRG")
  f2 <- fwb(d, lm_stat, R = 8L, verbose = FALSE, cluster = outside)

  expect_same_t(f2, f)

  #A non-factor cluster vector is accepted and converted.
  d[["gi"]] <- as.integer(d[["g"]])

  set.seed(78, "L'Ecuyer-CMRG")
  f3 <- fwb(d, lm_stat, R = 8L, verbose = FALSE, cluster = gi)

  expect_same_t(f3, f)
})

test_that("clusters and strata combine, with nesting enforced", {
  d <- fwb_test_df()

  set.seed(79, "L'Ecuyer-CMRG")
  expect_no_error({
    f <- fwb(d, w_identity_stat, R = 12L, verbose = FALSE, cluster = g,
             strata = s, simple = FALSE)
  })

  expect_output(print(f), "STRATIFIED FRACTIONAL WEIGHTED CLUSTER BOOTSTRAP")

  #The stratum totals still hold, and the weights are still constant within cluster.
  sums <- stratum_sums(f[["t"]], d[["s"]])

  expect_equal(as.vector(sums[1L, ]), as.numeric(table(d[["s"]])),
               tolerance = fwb_eps())

  expect_true(all(apply(f[["t"]], 1L, function(w) {
    all(tapply(w, d[["g"]], function(z) fwb:::all_the_same(z)))
  })))

  #A cluster that straddles two strata is rejected: the stratum a cluster belongs to
  #has to be well defined before a single weight can be drawn for it. The period-5
  #cycle below is deliberately coprime to `g`'s period of 12, so each cluster spans
  #several strata.
  d[["bad"]] <- factor(rep(c("a", "b", "c", "d", "e"), length.out = nrow(d)))

  expect_err(fwb(d, lm_stat, R = 8L, verbose = FALSE, cluster = g, strata = bad,
                 simple = FALSE),
             "must be completely nested within strata")
})

test_that("the `print()` header reflects strata and clusters", {
  d <- fwb_test_df()

  set.seed(80, "L'Ecuyer-CMRG")
  expect_output(print(fwb(d, lm_stat, R = 8L, verbose = FALSE)),
                "^FRACTIONAL WEIGHTED BOOTSTRAP")

  set.seed(80, "L'Ecuyer-CMRG")
  expect_output(print(fwb(d, lm_stat, R = 8L, verbose = FALSE, cluster = g)),
                "^FRACTIONAL WEIGHTED CLUSTER BOOTSTRAP")

  set.seed(80, "L'Ecuyer-CMRG")
  expect_output(print(fwb(d, lm_stat, R = 8L, verbose = FALSE, strata = s,
                          simple = FALSE)),
                "^STRATIFIED FRACTIONAL WEIGHTED BOOTSTRAP")
})
