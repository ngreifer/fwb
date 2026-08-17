#`summary()`, `confint()`, `tidy()`, and the simultaneous-inference machinery.
#
#`test-p-values.R` already checks the central property of the p-values -- that
#inverting a confidence interval at level `1 - p` reproduces the interval that put
#the null on its boundary. This file covers the surrounding structure: which columns
#appear for which arguments, that `summary()` and `confint()` and `fwb.ci()` agree,
#and the sup-t bands.

boot_for_summary <- function(n = 40L, R = 200L, seed = 91L) {
  d <- fwb_test_df(n)

  set.seed(seed, "L'Ecuyer-CMRG")
  fwb(d, lm_stat, R = R, verbose = FALSE, simple = FALSE)
}

test_that("`summary()` reports estimates, standard errors, and intervals", {
  f <- boot_for_summary()

  s <- summary(f)

  expect_s3_class(s, "summary.fwb")
  expect_true(is.matrix(s))
  expect_identical(rownames(s), names(f[["t0"]]))
  expect_identical(colnames(s), c("Estimate", "Std. Error", "CI 2.5 %", "CI 97.5 %"))

  expect_equal(s[, "Estimate"], f[["t0"]], tolerance = fwb_eps())
  expect_equal(s[, "Std. Error"], apply(f[["t"]], 2L, sd),
               tolerance = fwb_eps(), ignore_attr = TRUE)

  expect_identical(attr(s, "conf", TRUE), .95)
  expect_identical(attr(s, "ci.type", TRUE), "bc")
  expect_false(attr(s, "simultaneous", TRUE))
  expect_null(attr(s, "null", TRUE))

  #The level appears in the column labels.
  s90 <- summary(f, conf = .9)

  expect_identical(colnames(s90)[3:4], c("CI 5 %", "CI 95 %"))
})

test_that("`conf = 0` suppresses the interval columns", {
  f <- boot_for_summary()

  s <- summary(f, conf = 0)

  expect_identical(colnames(s), c("Estimate", "Std. Error"))
  expect_identical(attr(s, "conf", TRUE), 0)

  #With `p.value = TRUE` the p-value column is still produced.
  s <- summary(f, conf = 0, p.value = TRUE)

  expect_identical(colnames(s), c("Estimate", "Std. Error", "Pr(>|z|)"))

  #`ci.type = "wald"` adds a z statistic, before the p-value.
  s <- summary(f, conf = 0, p.value = TRUE, ci.type = "wald")

  expect_identical(colnames(s),
                   c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))

  expect_equal(s[, "z value"], (f[["t0"]] - 0) / apply(f[["t"]], 2L, sd),
               tolerance = fwb_eps(), ignore_attr = TRUE)

  #And with both intervals and a z statistic, the column order is
  #estimate, SE, interval, statistic, p-value -- which `print.summary.fwb()` relies
  #on to pick out `tst.ind`.
  s <- summary(f, p.value = TRUE, ci.type = "wald")

  expect_identical(colnames(s),
                   c("Estimate", "Std. Error", "CI 2.5 %", "CI 97.5 %",
                     "z value", "Pr(>|z|)"))
})

test_that("`null` shifts the hypothesis being tested", {
  f <- boot_for_summary()

  s0 <- summary(f, conf = 0, p.value = TRUE, ci.type = "wald", index = 2L)
  s1 <- summary(f, conf = 0, p.value = TRUE, ci.type = "wald", index = 2L,
                null = f[["t0"]][[2L]])

  expect_identical(attr(s1, "null", TRUE), f[["t0"]][[2L]])

  #Testing against the estimate itself gives z = 0 and p = 1.
  expect_equal(s1[1L, "z value"], 0, tolerance = fwb_eps(), ignore_attr = TRUE)
  expect_equal(s1[1L, "Pr(>|z|)"], 1, tolerance = fwb_eps(), ignore_attr = TRUE)

  expect_gt(s1[1L, "Pr(>|z|)"], s0[1L, "Pr(>|z|)"])
})

test_that("`index` selects statistics by name or position", {
  f <- boot_for_summary()

  by_pos <- summary(f, index = c(1L, 3L))
  by_name <- summary(f, index = names(f[["t0"]])[c(1L, 3L)])

  expect_identical(rownames(by_pos), names(f[["t0"]])[c(1L, 3L)])
  expect_equal(unname(by_pos), unname(by_name), tolerance = fwb_eps())

  #Duplicates are collapsed rather than repeated.
  expect_identical(nrow(summary(f, index = c(2L, 2L))), 1L)

  #`confint()` takes the same values through `parm`.
  expect_identical(rownames(confint(f, parm = 2L)), names(f[["t0"]])[2L])
  expect_identical(rownames(confint(f, parm = "x1")), "x1")

  #Missing `parm` means every statistic.
  expect_identical(nrow(confint(f)), length(f[["t0"]]))

  #The messages name the argument the *caller* used, which is `parm` in `confint()` and
  #`index` everywhere else. They used to report the offending value instead -- and after
  #`match()` had overwritten it, so a bad name was reported as "`NA_integer_`".
  expect_err(summary(f, index = 99L), "`index` must be between 1 and")
  expect_err(summary(f, index = "nope"), "all entries in `index` must be the names")

  expect_err(confint(f, parm = 99L), "`parm` must be between 1 and")
  expect_err(confint(f, parm = "nope"), "all entries in `parm` must be the names")
})

test_that("`summary()`, `confint()`, and `fwb.ci()` agree on every interval type", {
  f <- boot_for_summary()

  index <- 2L

  for (ty in c("wald", "norm", "cheap", "basic", "perc", "bc", "bca")) {
    ci <- suppressWarnings(get_ci(fwb.ci(f, conf = .9, type = ty,
                                         index = index))[[ty]])

    s <- suppressWarnings(summary(f, conf = .9, ci.type = ty, index = index))
    cf <- suppressWarnings(confint(f, parm = index, level = .9, ci.type = ty))

    expect_equal(unname(s[1L, 3:4]), unname(ci),
                 tolerance = fwb_eps(), info = ty)
    expect_equal(unname(cf[1L, ]), unname(ci),
                 tolerance = fwb_eps(), info = ty)
  }
})

test_that("`confint()` labels its columns with the requested level", {
  f <- boot_for_summary()

  cf <- confint(f, level = .8)

  expect_identical(colnames(cf), c("10 %", "90 %"))
  expect_identical(dim(cf), c(length(f[["t0"]]), 2L))
  expect_false(attr(cf, "simultaneous", TRUE))
})

test_that("`print.summary.fwb()` prints a coefficient table", {
  f <- boot_for_summary()

  out <- collapse_ws(utils::capture.output(print(summary(f, p.value = TRUE,
                                                         ci.type = "wald"))))

  expect_match(out, "Estimate", fixed = TRUE)
  expect_match(out, "Std. Error", fixed = TRUE)
  expect_match(out, "z value", fixed = TRUE)
  expect_match(out, "Signif. codes", fixed = TRUE)

  #No p-values means no significance codes.
  out <- collapse_ws(utils::capture.output(print(summary(f))))

  expect_no_match(out, "Signif. codes", fixed = TRUE)

  #`print()` returns its argument invisibly.
  s <- summary(f)

  utils::capture.output({
    res <- print(s)
  })

  expect_identical(res, s)
})

test_that("`tidy()` renames the columns to broom conventions", {
  f <- boot_for_summary()

  #Interval and p-value present.
  td <- generics::tidy(summary(f, p.value = TRUE, ci.type = "wald"))

  expect_s3_class(td, "tbl_df")
  expect_identical(names(td),
                   c("term", "estimate", "std.error", "conf.low", "conf.high",
                     "statistic", "p.value"))
  expect_identical(td[["term"]], names(f[["t0"]]))

  #The statistic names move into `term`, so the row names are the default sequence
  #rather than a duplicate of that column.
  expect_identical(rownames(td), as.character(seq_len(nrow(td))))

  #Interval only.
  td <- generics::tidy(summary(f))

  expect_identical(names(td),
                   c("term", "estimate", "std.error", "conf.low", "conf.high"))

  #p-value only, no z statistic (non-Wald interval type).
  td <- generics::tidy(summary(f, conf = 0, p.value = TRUE, ci.type = "perc"))

  expect_identical(names(td), c("term", "estimate", "std.error", "p.value"))

  #z statistic without an interval.
  td <- generics::tidy(summary(f, conf = 0, p.value = TRUE, ci.type = "wald"))

  expect_identical(names(td),
                   c("term", "estimate", "std.error", "statistic", "p.value"))
})

test_that("`summary()` validates its arguments", {
  f <- boot_for_summary()

  expect_err(summary(f, conf = 1), "greater than or equal to 0 and less than 1")
  expect_err(summary(f, conf = -.1), "greater than or equal to 0 and less than 1")
  expect_err(summary(f, ci.type = "nope"), "should be one of")
  expect_err(summary(f, ci.type = c("perc", "bc")), "must be a string")
  expect_err(summary(f, p.value = "yes"), "must be")
  expect_err(summary(f, p.value = TRUE, null = c(0, 1)), "must be a single number")

  #`level = 0` is accepted by `summary()` (as `conf`) but rejected by `confint()`,
  #even though the shared documentation says both may be set to 0 to suppress the
  #interval. `summary(conf = 0)` is the documented way to do it.
  expect_no_error({
    s <- summary(f, conf = 0)
  })

  expect_err(confint(f, level = 0), "must be between 0 and 1")
})

test_that("simultaneous Wald bands are wider than pointwise and stay ordered", {
  skip_if_not_installed("mvtnorm")

  f <- boot_for_summary()

  index <- seq_len(3L)

  pointwise <- confint(f, parm = index, level = .95, ci.type = "wald")
  simul <- confint(f, parm = index, level = .95, ci.type = "wald",
                   simultaneous = TRUE)

  expect_true(attr(simul, "simultaneous", TRUE))

  #A sup-t band must cover the whole vector, so each of its intervals is at least as
  #wide as the pointwise one -- and strictly wider once the estimates are not
  #perfectly correlated.
  expect_true(all(simul[, 1L] <= pointwise[, 1L]))
  expect_true(all(simul[, 2L] >= pointwise[, 2L]))
  expect_true(any(simul[, 2L] - simul[, 1L] > pointwise[, 2L] - pointwise[, 1L]))

  #And less conservative than Bonferroni, which is the claim in `?summary.fwb`.
  bonferroni <- confint(f, parm = index, level = 1 - .05 / length(index),
                        ci.type = "wald")

  expect_true(all(simul[, 1L] >= bonferroni[, 1L]))
  expect_true(all(simul[, 2L] <= bonferroni[, 2L]))

  #Simultaneous p-values are correspondingly larger.
  p_point <- summary(f, conf = 0, p.value = TRUE, ci.type = "wald",
                     index = index)[, "Pr(>|z|)"]
  p_simul <- summary(f, conf = 0, p.value = TRUE, ci.type = "wald",
                     index = index, simultaneous = TRUE)[, "Pr(>|z|)"]

  expect_true(all(p_simul >= p_point - fwb_eps()))
})

test_that("simultaneous percentile bands are wider than pointwise", {
  f <- boot_for_summary()

  index <- seq_len(3L)

  pointwise <- confint(f, parm = index, level = .9, ci.type = "perc")
  simul <- suppressWarnings(confint(f, parm = index, level = .9,
                                    ci.type = "perc", simultaneous = TRUE))

  expect_true(attr(simul, "simultaneous", TRUE))
  expect_true(all(simul[, 1L] <= pointwise[, 1L] + fwb_eps()))
  expect_true(all(simul[, 2L] >= pointwise[, 2L] - fwb_eps()))
})

test_that("the simultaneous percentile level attains its nominal coverage", {
  #The adjusted level is the narrowest pointwise level whose intervals jointly contain at
  #least `level` of the bootstrap draws. It is read off the order statistics rather than
  #searched for, so every assertion here is checked against `simultaneous_coverage()` --
  #which builds the intervals and counts, i.e., the definition rather than the shortcut.
  f <- boot_for_summary(R = 300L)

  R <- nrow(f[["t"]])

  for (level in c(.8, .9, .95)) {
    for (k in 2:3) {
      index <- seq_len(k)
      info <- sprintf("level = %s, k = %s", level, k)

      l <- fwb:::simultaneous_ci_level(f, level, index, "perc")

      #The property that matters: the band attains at least what was asked for.
      expect_gte(simultaneous_coverage(f, index, l), level, label = info)

      #And it is the narrowest level that does: one whole order-statistic step below it,
      #the target is missed. This is what rules out simply returning the conservative
      #Bonferroni level, which would also attain the target.
      expect_lt(simultaneous_coverage(f, index, l - 2 / (R + 1)), level)
    }
  }
})

test_that("the simultaneous percentile level matches an exhaustive search", {
  #An independent check on the algebra behind the one-pass computation: the achievable
  #levels are `|2i/(R + 1) - 1|` for integer `i`, so the answer can be found by evaluating
  #the coverage at every one of them and taking the smallest that attains the target. Slow,
  #but it shares nothing with the implementation except the definition of coverage.
  f <- boot_for_summary(R = 120L)

  R <- nrow(f[["t"]])

  candidates <- sort(unique(abs(2 * seq_len(R) / (R + 1) - 1)))

  for (level in c(.85, .95)) {
    for (k in 2:3) {
      index <- seq_len(k)

      covered <- vapply(candidates,
                        function(l) simultaneous_coverage(f, index, l) >= level,
                        logical(1L))

      brute <- candidates[which(covered)[1L]]

      expect_equal(fwb:::simultaneous_ci_level(f, level, index, "perc"), brute,
                   tolerance = 1e-8,
                   label = sprintf("level = %s, k = %s", level, k))
    }
  }
})

test_that("the exact level is returned a few ulps wide, on purpose", {
  #At the exact level the binding draw sits *on* an endpoint. `compute_ci()` is not handed
  #`p`; it is handed the level and recomputes `p = (1 - level) / 2`, and because the level
  #is near 1 that round trip loses about an ulp of `p` -- enough to flip
  #`floor((R + 1) * p)` and drop the binding draw, leaving coverage one draw short.
  #
  #So the returned level is nudged a few ulps wider. The test is that this is invisible in
  #the interval but decisive for the coverage.
  f <- boot_for_summary(R = 300L)

  index <- 1:3
  level <- .9

  l <- fwb:::simultaneous_ci_level(f, level, index, "perc")
  exact <- l - 16 * .Machine$double.eps

  expect_gte(simultaneous_coverage(f, index, l), level)

  limits <- function(conf) {
    m <- suppressWarnings(fwb:::compute_ci("perc", t = f[["t"]], t0 = f[["t0"]],
                                           conf = conf, index = index, boot.out = f))

    m[, (ncol(m) - 1L):ncol(m), drop = FALSE]
  }

  #The endpoints move by far less than the data's own precision.
  expect_lt(max(abs(limits(l) - limits(exact))), 1e-10)
})

test_that("simultaneous percentile p-values invert the bands", {
  #`simultaneous_ci_level()` inverts the coverage function and `simultaneous_p_value()`
  #evaluates it; both read it off the same entry levels, which is what keeps them
  #consistent. That shortcut has to equal the interval-by-interval computation, or the
  #p-values would no longer invert the intervals.
  f <- boot_for_summary(R = 300L)

  index <- 1:3

  entry <- fwb:::sup_t_entry_levels(f[["t"]][, index, drop = FALSE])

  levels <- seq(.05, .99, by = .01)

  expect_equal(vapply(levels, function(l) mean(entry <= l), numeric(1L)),
               vapply(levels, function(l) simultaneous_coverage(f, index, l), numeric(1L)),
               tolerance = 0)

  #And the p-value reported for a null on the band's own boundary is one minus the
  #coverage the band attains.
  l <- fwb:::simultaneous_ci_level(f, .9, index, "perc")

  p <- fwb:::simultaneous_p_value(f, rep(1 - l, length(index)), index, "perc")

  expect_equal(p, rep(1 - simultaneous_coverage(f, index, l), length(index)),
               tolerance = fwb_eps())
})

test_that("simultaneous inference rejects non-finite estimates", {
  #This used to produce `NA` without saying anything: the coverage comparison ran against
  #`NA` endpoints, and `cov()` on the Wald side returned `NA`.
  f <- boot_for_summary(R = 300L)

  is.na(f[["t"]][3L, 2L]) <- TRUE

  expect_err(confint(f, parm = 1:3, ci.type = "perc", simultaneous = TRUE),
             "non-finite")
  expect_err(summary(f, conf = 0, p.value = TRUE, ci.type = "perc",
                     simultaneous = TRUE),
             "non-finite")
})

test_that("simultaneous inference is restricted to `\"wald\"` and `\"perc\"`", {
  f <- boot_for_summary()

  for (ty in c("norm", "cheap", "basic", "bc", "bca")) {
    expect_err(summary(f, ci.type = ty, simultaneous = TRUE),
               "simultaneous inference can only be used")
    expect_err(confint(f, level = .95, ci.type = ty, simultaneous = TRUE),
               "simultaneous inference can only be used")
  }

  #Simultaneous Wald needs a level above .5, because the sup-t critical value is
  #derived from a two-sided multivariate normal quantile.
  expect_err(summary(f, conf = .4, ci.type = "wald", simultaneous = TRUE),
             "must be greater than .5")
  expect_err(confint(f, level = .4, ci.type = "wald", simultaneous = TRUE),
             "must be greater than .5")

  #With a single statistic there is nothing to adjust, so `simultaneous` is silently
  #ignored rather than being an error.
  s <- summary(f, index = 1L, simultaneous = TRUE)

  expect_false(attr(s, "simultaneous", TRUE))

  cf <- confint(f, parm = 1L, simultaneous = TRUE)

  expect_false(attr(cf, "simultaneous", TRUE))
})

test_that("simultaneous Wald inference handles a zero-variance estimate", {
  #`simultaneous_p_value()` drops zero-variance statistics from the correlation
  #matrix, so the `vapply()` that computes the p-values has to be over the *subset* of
  #`z` -- computing it over the full `z` and assigning into `p[!zeros]` recycles, and
  #the surviving p-values land in the wrong slots.
  skip_if_not_installed("mvtnorm")

  d <- fwb_test_df()

  const_stat <- function(data, w) {
    c(coef(lm(y ~ x1, data = data, weights = w)), const = 1)
  }

  set.seed(92, "L'Ecuyer-CMRG")
  f <- fwb(d, const_stat, R = 200L, verbose = FALSE, simple = FALSE)

  expect_no_warning({
    s <- summary(f, conf = 0, p.value = TRUE, ci.type = "wald",
                 simultaneous = TRUE)
  })

  #The two real coefficients keep sensible p-values; only the constant is degenerate.
  p_point <- summary(f, conf = 0, p.value = TRUE, ci.type = "wald")[, "Pr(>|z|)"]

  expect_true(all(s[1:2, "Pr(>|z|)"] >= p_point[1:2] - fwb_eps()))
})

test_that("`coef()` and `vcov()` methods read `t0` and `t`", {
  f <- boot_for_summary()

  expect_identical(coef(f), f[["t0"]])
  expect_equal(vcov(f), cov(f[["t"]]), tolerance = fwb_eps())

  #The standard errors in `summary()` are the square roots of the `vcov()` diagonal.
  expect_equal(summary(f)[, "Std. Error"], sqrt(diag(vcov(f))),
               tolerance = fwb_eps(), ignore_attr = TRUE)
})
