#Confidence intervals: `fwb.ci()`, `get_ci()`, and the interval arithmetic in
#`compute_ci()`.
#
#`test-boot.R` already checks agreement with `boot::boot.ci()` for the four interval
#types the two packages share. What is left uncovered, and is what this file adds, is
#the *fwb*-only machinery: the `"wald"`, `"bc"`, and `"cheap"` intervals, the
#identities that tie the types together, and the argument handling.

boot_for_ci <- function(n = 40L, R = 200L, seed = 81L) {
  d <- fwb_test_df(n)

  set.seed(seed, "L'Ecuyer-CMRG")
  fwb(d, lm_stat, R = R, verbose = FALSE, simple = FALSE)
}

test_that("`fwb.ci()` returns one component per requested type", {
  f <- boot_for_ci()

  types <- c("perc", "bc", "wald", "norm", "cheap", "basic", "bca")

  ci <- fwb.ci(f, type = types, index = 2L)

  expect_s3_class(ci, c("fwbci", "bootci"))
  expect_true(all(types %in% names(ci)))
  expect_identical(ci[["R"]], f[["R"]])
  expect_equal(ci[["t0"]], f[["t0"]][[2L]], tolerance = fwb_eps(),
               ignore_attr = TRUE)
  expect_identical(attr(ci, "conf", TRUE), .95)

  #`"wald"`, `"norm"`, and `"cheap"` report level and the two limits; the
  #order-statistic types additionally report the two ranks used.
  for (ty in c("wald", "norm", "cheap")) {
    expect_identical(dim(ci[[ty]]), c(1L, 3L), info = ty)
  }

  for (ty in c("perc", "bc", "basic", "bca")) {
    expect_identical(dim(ci[[ty]]), c(1L, 5L), info = ty)
  }

  #`type = "all"` is the full set when there are no clusters.
  ci_all <- fwb.ci(f, type = "all", index = 2L)

  expect_true(all(types %in% names(ci_all)))
})

test_that("`type = \"all\"` drops BCa when clusters are present", {
  d <- fwb_test_df()

  set.seed(82, "L'Ecuyer-CMRG")
  f <- fwb(d, lm_stat, R = 200L, verbose = FALSE, simple = FALSE, cluster = g)

  ci <- fwb.ci(f, type = "all", index = 2L)

  expect_false("bca" %in% names(ci))
  expect_true(all(c("perc", "bc", "wald", "norm", "cheap", "basic") %in% names(ci)))
})

test_that("the interval types satisfy their defining identities", {
  f <- boot_for_ci()

  index <- 2L
  t <- f[["t"]][, index]
  t0 <- f[["t0"]][[index]]
  conf <- .95
  z <- qnorm(1 - (1 - conf) / 2)

  ci <- fwb.ci(f, conf = conf, type = c("wald", "norm", "basic", "perc", "bc"),
               index = index)

  #Wald: t0 +/- z * sd, with no bias correction at all.
  expect_equal(as.vector(ci[["wald"]][1L, 2:3]),
               t0 + c(-1, 1) * z * sd(t),
               tolerance = fwb_eps())

  #Normal: the same, after subtracting the bootstrap bias.
  expect_equal(as.vector(ci[["norm"]][1L, 2:3]),
               t0 - (mean(t) - t0) + c(-1, 1) * z * sd(t),
               tolerance = fwb_eps())

  #Basic is the percentile interval reflected about t0, so the two must sum to 2*t0
  #limit by limit (in reverse order).
  expect_equal(as.vector(ci[["basic"]][1L, 4:5]),
               2 * t0 - rev(as.vector(ci[["perc"]][1L, 4:5])),
               tolerance = fwb_eps())

  #`"bc"` shifts the percentile ranks by the median bias. When t0 sits at the median
  #of the bootstrap distribution the shift is zero and the two coincide -- which is
  #what `?fwb.ci` claims.
  bias_frac <- mean(t < t0)

  if (abs(bias_frac - 0.5) < 1e-12) {
    expect_equal(as.vector(ci[["bc"]][1L, 4:5]),
                 as.vector(ci[["perc"]][1L, 4:5]),
                 tolerance = fwb_eps())
  }

  #Constructed exactly: recentering the replicates on their own median makes the
  #median-bias correction vanish.
  f2 <- f
  f2[["t0"]][index] <- median(t)

  ci2 <- fwb.ci(f2, conf = conf, type = c("perc", "bc"), index = index)

  expect_equal(as.vector(ci2[["bc"]][1L, 4:5]),
               as.vector(ci2[["perc"]][1L, 4:5]),
               tolerance = 1e-6)
})

test_that("the `\"cheap\"` interval matches Lam's formula and works at tiny `R`", {
  f <- boot_for_ci()

  index <- 2L
  t <- f[["t"]][, index]
  t0 <- f[["t0"]][[index]]
  conf <- .9

  ci <- fwb.ci(f, conf = conf, type = "cheap", index = index)

  s_star <- sqrt(mean((t - t0)^2))
  tcrit <- qt(1 - (1 - conf) / 2, df = length(t))

  expect_equal(as.vector(ci[["cheap"]][1L, 2:3]),
               t0 + c(-1, 1) * s_star * tcrit,
               tolerance = fwb_eps())

  #The point of the cheap interval is that it stays defined with almost no
  #replicates -- the reference calls out R = 1. Every other type needs order
  #statistics or a bootstrap SD and cannot.
  d <- fwb_test_df()

  set.seed(83, "L'Ecuyer-CMRG")
  f1 <- fwb(d, lm_stat, R = 1L, verbose = FALSE)

  for (route in list(function(x) confint(x, level = conf, ci.type = "cheap"),
                     function(x) {
                       m <- fwb.ci(x, conf = conf, type = "cheap")[["cheap"]]
                       m[, 2:3, drop = FALSE]
                     })) {
    ci1 <- route(f1)

    expect_true(all(is.finite(ci1)))
    expect_true(all(ci1[, 1L] < ci1[, 2L]))
  }

  #And the `summary()` route reports the p-value against a t distribution rather
  #than a normal.
  s <- summary(f1, ci.type = "cheap", p.value = TRUE)

  expect_true("Pr(>|t|)" %in% colnames(s))
  expect_false("Pr(>|z|)" %in% colnames(s))

  #With one replicate the SD is undefined, so `Std. Error` is `NA` while the
  #interval is not.
  expect_true(all(is.na(s[, "Std. Error"])))
  expect_true(all(is.finite(s[, 3:4])))
})

test_that("intervals widen with the confidence level and bracket the estimate", {
  f <- boot_for_ci()

  index <- 2L
  t0 <- f[["t0"]][[index]]

  for (ty in c("wald", "norm", "cheap", "basic", "perc", "bc", "bca")) {
    narrow <- suppressWarnings(get_ci(fwb.ci(f, conf = .80, type = ty,
                                             index = index))[[ty]])
    wide <- suppressWarnings(get_ci(fwb.ci(f, conf = .99, type = ty,
                                           index = index))[[ty]])

    expect_lt(narrow[["L"]], narrow[["U"]])
    expect_lt(wide[["L"]], wide[["U"]])

    expect_lt(wide[["L"]], narrow[["L"]])
    expect_gt(wide[["U"]], narrow[["U"]])

    #All seven types must contain the point estimate for a statistic this
    #well-behaved; an interval that did not would signal a sign or ordering error.
    expect_lte(narrow[["L"]], t0)
    expect_gte(narrow[["U"]], t0)
  }
})

test_that("`hinv` transforms the limits and leaves the ranks alone", {
  f <- boot_for_ci()

  index <- 2L

  plain <- fwb.ci(f, type = c("perc", "norm"), index = index)
  transformed <- fwb.ci(f, type = c("perc", "norm"), index = index, hinv = exp)

  expect_equal(exp(as.vector(plain[["perc"]][1L, 4:5])),
               as.vector(transformed[["perc"]][1L, 4:5]),
               tolerance = fwb_eps())

  #The order statistics used are a property of the replicates, not of the scale.
  expect_equal(as.vector(plain[["perc"]][1L, 2:3]),
               as.vector(transformed[["perc"]][1L, 2:3]),
               tolerance = fwb_eps())

  expect_equal(exp(as.vector(plain[["norm"]][1L, 2:3])),
               as.vector(transformed[["norm"]][1L, 2:3]),
               tolerance = fwb_eps())

  #`t0` is reported on the transformed scale too.
  expect_equal(transformed[["t0"]], exp(f[["t0"]][[index]]),
               tolerance = fwb_eps(), ignore_attr = TRUE)
})

test_that("`get_ci()` extracts limits from `fwb.ci()` and `boot.ci()` alike", {
  f <- boot_for_ci()

  ci <- fwb.ci(f, conf = .9, type = c("perc", "bc", "wald"), index = 2L)

  g <- get_ci(ci)

  expect_named(g, c("perc", "bc", "wald"))
  expect_identical(attr(g, "conf", TRUE), .9)

  for (ty in names(g)) {
    expect_named(g[[ty]], c("L", "U"))
    expect_equal(unname(g[[ty]]),
                 unname(fwb:::.tail(ci[[ty]][1L, ], 2L)),
                 tolerance = fwb_eps(), info = ty)
  }

  #A subset can be requested by name.
  expect_named(get_ci(ci, "bc"), "bc")
  expect_named(get_ci(ci, c("wald", "perc")), c("wald", "perc"))

  expect_err(get_ci(ci, "bca"), "should be at least one of")
  expect_err(get_ci(f), "must inherit from class")

  #*boot* names two of its interval types differently; `get_ci()` renames them so the
  #two packages' output can be compared directly.
  skip_if_not_installed("boot")

  d <- fwb_test_df(40L)

  set.seed(84)
  b <- boot::boot(d, function(dat, i) coef(lm(y ~ x1 + x2, data = dat[i, ])),
                  R = 200L)

  bci <- boot::boot.ci(b, type = c("norm", "perc"), index = 2L)

  gb <- get_ci(bci)

  #*boot* returns its components in its own order, which is not the requested one.
  expect_setequal(names(gb), c("norm", "perc"))
  expect_named(gb[["norm"]], c("L", "U"))
})

test_that("`fwb.ci()` refuses BCa when it cannot be computed", {
  d <- fwb_test_df(60L)

  #Fewer replicates than units: the influence regression is underdetermined.
  set.seed(85, "L'Ecuyer-CMRG")
  few <- fwb(d, lm_stat, R = 40L, verbose = FALSE, simple = FALSE)

  expect_err(fwb.ci(few, type = "bca"),
             "fewer bootstrap replications than units")

  #Requested alongside other types, BCa is dropped with a warning rather than
  #erroring, so the rest of the output survives.
  expect_wrn(ci <- fwb.ci(few, type = c("perc", "bca")),
             "fewer bootstrap replications than units")

  expect_true("perc" %in% names(ci))
  expect_true(fwb:::is_null(ci[["bca"]]))

  #Clusters.
  set.seed(85, "L'Ecuyer-CMRG")
  clus <- fwb(d, lm_stat, R = 200L, verbose = FALSE, simple = FALSE, cluster = g)

  expect_err(fwb.ci(clus, type = "bca"), "cannot be used with clusters")
})

test_that("BCa intervals work with a random `statistic` and `simple = TRUE`", {
  #This combination used to be refused, because recovering the weights meant replaying a
  #stream that the statistic's own draws had interleaved with. Each replicate's weights
  #now come from a stream recorded before the run, so the influence values can always be
  #computed and the refusal is gone -- along with the malformed cli message that used to
  #report it as `Error: Expecting '}'`.
  d <- fwb_test_df(40L)

  set.seed(86, "L'Ecuyer-CMRG")
  f <- fwb(d, random_lm_stat, R = 200L, verbose = FALSE, simple = TRUE)

  expect_no_error({
    ci <- fwb.ci(f, type = "bca", index = 2L)
  })

  expect_s3_class(ci, "fwbci")
  expect_true(all(is.finite(ci[["bca"]][1L, 4:5])))

  #The remaining BCa restrictions still hold and still report cleanly.
  set.seed(86, "L'Ecuyer-CMRG")
  clus <- fwb(d, random_lm_stat, R = 200L, verbose = FALSE, cluster = g)

  expect_err(fwb.ci(clus, type = "bca"), "cannot be used with clusters")
})

test_that("`fwb.ci()` validates its arguments", {
  f <- boot_for_ci(R = 200L)

  expect_err(fwb.ci(f, conf = 1.5), "must be between 0 and 1")
  expect_err(fwb.ci(f, conf = 0), "must be between 0 and 1")
  expect_err(fwb.ci(f, conf = c(.9, .95)), "must be a single number")
  expect_err(fwb.ci(f, type = "nope"), "should be at least one of")
  expect_err(fwb.ci(f, index = 1:2), "`index` must have length one")
  expect_err(fwb.ci(f, index = 99), "`index` must be between 1 and")
  expect_err(fwb.ci(f, index = "nope"), "all entries in `index` must be the names")
  expect_err(fwb.ci(fwb_test_df()), "must inherit from class")

  #`h` is accepted for `boot.ci()` compatibility but only the identity is supported.
  expect_err(fwb.ci(f, h = exp), "can only be")
  expect_no_error({
    ci <- fwb.ci(f, h = identity)
  })

  #`NA` and non-finite replicates make every interval undefined.
  fna <- f
  is.na(fna[["t"]][3L, 1L]) <- TRUE

  expect_err(fwb.ci(fna, index = 1L), "cannot calculate confidence intervals")

  finf <- f
  finf[["t"]][3L, 1L] <- Inf

  expect_err(fwb.ci(finf, index = 1L), "non-finite")

  #A degenerate statistic.
  d <- fwb_test_df()

  set.seed(87, "L'Ecuyer-CMRG")
  fconst <- fwb(d, function(data, w) c(a = 1), R = 30L, verbose = FALSE)

  expect_err(fwb.ci(fconst), "all estimates are equal")
})

test_that("extreme confidence levels warn about the order statistics used", {
  f <- boot_for_ci(R = 200L)

  for (ty in c("perc", "bc", "basic")) {
    expect_wrn(ci <- fwb.ci(f, conf = .9999, type = ty, index = 2L),
               "extreme order statistics used as endpoints")
  }

  #The closed-form types have no order statistics to run out of.
  for (ty in c("wald", "norm", "cheap")) {
    expect_no_warning({
      ci <- fwb.ci(f, conf = .9999, type = ty, index = 2L)
    })
  }
})

test_that("`print.fwbci()` reports the types, level, and scale", {
  f <- boot_for_ci(R = 200L)

  ci <- fwb.ci(f, conf = .9, type = c("perc", "bc", "wald", "basic"), index = 2L)

  out <- utils::capture.output(print(ci))
  txt <- collapse_ws(out)

  expect_match(txt, "BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS", fixed = TRUE)
  expect_match(txt, "Based on 200 bootstrap replicates", fixed = TRUE)
  expect_match(txt, "90%", fixed = TRUE)
  expect_match(txt, "Percentile", fixed = TRUE)
  expect_match(txt, "BC Percentile", fixed = TRUE)
  expect_match(txt, "Wald", fixed = TRUE)
  expect_match(txt, "Basic", fixed = TRUE)
  expect_match(txt, "Calculations and Intervals on Original Scale", fixed = TRUE)

  #`hinv` in the original call is reported as such.
  ci2 <- fwb.ci(f, type = "perc", index = 2L, hinv = exp)

  expect_match(collapse_ws(utils::capture.output(print(ci2))),
               "Intervals Transformed", fixed = TRUE)

  #`print()` returns its input invisibly.
  utils::capture.output({
    res <- print(ci)
  })

  expect_identical(res, ci)
})

test_that("`norm_inter()` interpolates like `boot`'s version", {
  #`norm_inter()` uses a partial sort for speed, so it only guarantees that the two
  #order statistics it needs are in place. If the partial-sort indices ever fell out
  #of step with the indices used for interpolation, the result would be a plausible
  #but wrong quantile -- silently.
  set.seed(88, "L'Ecuyer-CMRG")
  ti <- rnorm(199L)

  alpha <- c(.025, .1, .5, .9, .975)

  got <- fwb:::norm_inter(ti, alpha)

  expect_identical(dim(got), c(length(alpha), 2L))

  #Recompute without the partial sort.
  sorted <- sort(ti)
  R <- length(ti)
  k <- floor((R + 1) * alpha)

  z_k <- qnorm(k / (R + 1))
  z_k1 <- qnorm((k + 1) / (R + 1))
  wt <- (qnorm(alpha) - z_k) / (z_k1 - z_k)

  expect_equal(got[, 2L], (1 - wt) * sorted[k] + wt * sorted[k + 1L],
               tolerance = fwb_eps())

  #Already-sorted input must take the same path.
  expect_equal(fwb:::norm_inter(sorted, alpha), got, tolerance = fwb_eps())

  skip_if_not_installed("boot")

  expect_equal(unname(got), unname(boot:::norm.inter(ti, alpha)),
               tolerance = fwb_eps())
})

test_that("`norm_inter_invert()` inverts `norm_inter()`", {
  set.seed(89, "L'Ecuyer-CMRG")
  ti <- rnorm(199L)

  #Round-tripping a quantile through the interpolation and back must land where it
  #started. This is the identity the p-values in `summary()` rely on.
  for (alpha in c(.05, .25, .5, .75, .95)) {
    q <- fwb:::norm_inter(ti, alpha)[1L, 2L]

    expect_equal(fwb:::norm_inter_invert(ti, q), alpha, tolerance = 1e-6)
  }

  #Outside the support the inverse saturates rather than extrapolating.
  expect_identical(fwb:::norm_inter_invert(ti, min(ti) - 1), 0)
  expect_identical(fwb:::norm_inter_invert(ti, max(ti) + 1), 1)
})

test_that("interval types that need more than one replicate say so", {
  #Each type used to fail its own way below its minimum `R`, and none of them said why:
  #`"wald"`/`"norm"` returned `NA` limits silently, `"perc"`/`"basic"` errored from inside
  #`norm_inter()` with "subscript out of bounds", and `"bc"`/`"bca"` with "missing value
  #where TRUE/FALSE needed". The check lives in `compute_ci()`/`invert_ci()` rather than
  #`fwb.ci()` because `summary()` and `confint()` reach those directly.
  d <- fwb_test_df(40L)

  set.seed(90, "L'Ecuyer-CMRG")
  f1 <- fwb(d, lm_stat, R = 1L, verbose = FALSE)

  #`"cheap"` is the one that works at `R = 1` -- that is the point of it -- and it works
  #through all three interfaces.
  expect_no_error({
    ci <- fwb.ci(f1, type = "cheap")
  })

  expect_true(all(is.finite(confint(f1, ci.type = "cheap"))))
  expect_true(all(is.finite(summary(f1, ci.type = "cheap")[, 3:4])))
  expect_true(all(is.finite(summary(f1, conf = 0, p.value = TRUE,
                                    ci.type = "cheap")[, "Pr(>|t|)"])))

  #Every other type says what it needs, in each interface, for intervals and p-values
  #alike.
  for (ty in c("wald", "norm", "basic", "perc", "bc")) {
    expect_err(fwb.ci(f1, type = ty), "requires at least 2 bootstrap replications")
    expect_err(confint(f1, ci.type = ty), "requires at least 2 bootstrap replications")
    expect_err(summary(f1, ci.type = ty), "requires at least 2 bootstrap replications")
    expect_err(summary(f1, conf = 0, p.value = TRUE, ci.type = ty),
               "requires at least 2 bootstrap replications")
  }

  #The message names the type, so a `type = "all"` request is diagnosable.
  expect_err(fwb.ci(f1, type = "perc"), "\"perc\"")

  #Two replicates are enough for all of them except BCa, which needs more replications
  #than units.
  set.seed(90, "L'Ecuyer-CMRG")
  f2 <- fwb(d, lm_stat, R = 2L, verbose = FALSE)

  for (ty in c("cheap", "wald", "norm", "basic", "perc", "bc")) {
    expect_no_error({
      ci <- suppressWarnings(confint(f2, parm = 2L, ci.type = ty))
    })

    expect_true(all(is.finite(ci)), info = ty)
  }
})

test_that("`fwb.ci()` keeps the types it can compute and drops the rest", {
  #`type = "all"` at a small `R` should return what it can rather than failing on the
  #first type it cannot -- the same courtesy `fwb.ci()` already extended to BCa.
  d <- fwb_test_df(40L)

  set.seed(90, "L'Ecuyer-CMRG")
  f1 <- fwb(d, lm_stat, R = 1L, verbose = FALSE)

  expect_wrn(ci <- fwb.ci(f1, type = "all"),
             "more bootstrap replications are needed")

  expect_identical(intersect(names(ci), fwb:::.allowed_ci.types()), "cheap")

  #A partly computable request behaves the same way, and the warning names what went.
  expect_wrn(ci <- fwb.ci(f1, type = c("cheap", "perc", "wald")),
             "\"perc\"")

  expect_identical(intersect(names(ci), fwb:::.allowed_ci.types()), "cheap")

  #With nothing computable there is no output to return, so it errors instead.
  expect_err(fwb.ci(f1, type = c("perc", "wald")),
             "requires at least 2 bootstrap replications")

  #At a workable `R`, nothing is dropped and no warning is raised about it.
  set.seed(90, "L'Ecuyer-CMRG")
  f <- fwb(d, lm_stat, R = 200L, verbose = FALSE)

  expect_no_warning({
    ci <- fwb.ci(f, type = "all")
  }, message = "more bootstrap replications are needed")

  expect_setequal(intersect(names(ci), fwb:::.allowed_ci.types()),
                  fwb:::.allowed_ci.types())
})

test_that("`summary()` and `confint()` enforce BCa's requirements", {
  #`fwb.ci()` has always checked that BCa has more replications than units and no
  #clusters. `summary()` and `confint()` reach `compute_ci()` without passing through it,
  #so they used to compute an interval from an underdetermined influence regression and
  #report it without comment.
  d <- fwb_test_df(40L)

  set.seed(91, "L'Ecuyer-CMRG")
  few <- fwb(d, lm_stat, R = 30L, verbose = FALSE)

  expect_err(confint(few, parm = 2L, ci.type = "bca"),
             "more bootstrap replications than units")
  expect_err(summary(few, ci.type = "bca"),
             "more bootstrap replications than units")
  expect_err(summary(few, conf = 0, p.value = TRUE, ci.type = "bca"),
             "more bootstrap replications than units")

  set.seed(91, "L'Ecuyer-CMRG")
  clus <- fwb(d, lm_stat, R = 200L, verbose = FALSE, cluster = g)

  expect_err(confint(clus, parm = 2L, ci.type = "bca"), "cannot be used with clusters")
  expect_err(summary(clus, ci.type = "bca"), "cannot be used with clusters")

  #And with enough replications and no clusters it goes through.
  set.seed(91, "L'Ecuyer-CMRG")
  ok <- fwb(d, lm_stat, R = 200L, verbose = FALSE)

  expect_no_error({
    cf <- confint(ok, parm = 2L, ci.type = "bca")
  })

  expect_true(all(is.finite(cf)))
})
