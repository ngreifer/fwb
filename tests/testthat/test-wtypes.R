#Properties of the weight generators.
#
#`?fwb` makes specific distributional promises about each `wtype` -- mean, variance,
#and skewness -- and those promises are what justify the advice about which type to
#pick. They are also easy to break silently: a missing scale factor changes the
#variance without changing anything that would make a test of `fwb()` itself fail.

test_that("`.w_types()` and the documented set agree", {
  expect_identical(fwb:::.w_types(),
                   c("exp", "multinom", "poisson", "mammen", "beta", "power"))
})

test_that("each `wtype` has the documented mean, variance, and skewness", {
  skip_on_cran()

  #A large draw so the Monte Carlo error is small relative to the tolerances below,
  #but cheap: these are `rexp()`/`rbeta()`/`sample()` calls, not model fits.
  n <- 200L
  R <- 2000L

  #Skewness is the moment most sensitive to an error in the scaling constants, and
  #it is the reason `"mammen"`, `"beta"`, and `"power"` exist at all.
  expected <- list(exp      = c(mean = 1, var = 1, skew = 2),
                   multinom = c(mean = 1, var = 1, skew = 1),
                   poisson  = c(mean = 1, var = 1, skew = 1),
                   mammen   = c(mean = 1, var = 1, skew = 1),
                   beta     = c(mean = 1, var = 1, skew = 1),
                   power    = c(mean = 1, var = 1, skew = 2 * (sqrt(2) - 1)))

  for (wt in names(expected)) {
    gen <- fwb:::make_gen_weights(wt)

    expect_identical(attr(gen, "wtype", TRUE), wt)

    set.seed(50, "L'Ecuyer-CMRG")
    w <- as.vector(gen(n, R, NULL))

    got <- c(mean = mean(w),
             var = var(w),
             skew = mean((w - mean(w))^3) / sd(w)^3)

    expect_equal(got, expected[[wt]], tolerance = 0.05, info = wt)
  }
})

test_that("integer weight types produce non-negative integers", {
  for (wt in c("multinom", "poisson")) {
    gen <- fwb:::make_gen_weights(wt)

    set.seed(51, "L'Ecuyer-CMRG")
    w <- gen(40L, 20L, NULL)

    expect_true(all(w == round(w)), info = wt)
    expect_true(all(w >= 0), info = wt)

    #Both types are supposed to be able to zero a unit out -- that is what `drop0`
    #exists for -- so a generator that never did would make `drop0` untestable.
    expect_true(any(w == 0), info = wt)
  }

  #`"multinom"` is a tabulation of `n` draws with replacement, so each row must sum
  #to exactly `n`. That is what makes it identical to `boot::boot(stype = "f")`.
  gen <- fwb:::make_gen_weights("multinom")

  set.seed(52, "L'Ecuyer-CMRG")
  w <- gen(40L, 20L, NULL)

  expect_equal(rowSums(w), rep(40, 20L))
})

test_that("continuous weight types are strictly positive", {
  #No zeros means `drop0` has nothing to do, which is why `fwb()` forces it to
  #`FALSE` for these types.
  for (wt in c("exp", "mammen", "beta", "power")) {
    gen <- fwb:::make_gen_weights(wt)

    set.seed(53, "L'Ecuyer-CMRG")
    w <- gen(40L, 50L, NULL)

    expect_true(all(w > 0), info = wt)
  }

  #`"mammen"` takes exactly two values before normalization.
  gen <- fwb:::make_gen_weights("mammen")

  set.seed(54, "L'Ecuyer-CMRG")
  w <- gen(400L, 1L, NULL)

  expect_length(unique(round(as.vector(w), 10L)), 2L)
})

test_that("weights are normalized to mean 1 within each stratum", {
  n <- 90L
  strata <- factor(rep(c("a", "b"), times = c(30L, 60L)))

  #Every type except `"poisson"` rescales within stratum, so each stratum's weights
  #average 1 and the total weight equals the stratum size. This is what "essentially
  #a separate bootstrap within each stratum" means numerically.
  for (wt in c("exp", "multinom", "mammen", "beta", "power")) {
    gen <- fwb:::make_gen_weights(wt)

    set.seed(55, "L'Ecuyer-CMRG")
    w <- gen(n, 20L, strata)

    expect_equal(rowSums(w[, strata == "a", drop = FALSE]), rep(30, 20L),
                 tolerance = fwb_eps(), info = wt)
    expect_equal(rowSums(w[, strata == "b", drop = FALSE]), rep(60, 20L),
                 tolerance = fwb_eps(), info = wt)
  }
})

test_that("`\"poisson\"` ignores strata", {
  n <- 60L
  strata <- factor(rep(c("a", "b"), each = 30L))

  gen <- fwb:::make_gen_weights("poisson")

  set.seed(56, "L'Ecuyer-CMRG")
  with_strata <- gen(n, 10L, strata)

  set.seed(56, "L'Ecuyer-CMRG")
  without <- gen(n, 10L, NULL)

  expect_identical(with_strata, without)
})

test_that("a single-level stratum factor behaves like no strata", {
  n <- 40L
  strata <- factor(rep("only", n))

  for (wt in c("exp", "multinom", "mammen", "beta", "power")) {
    gen <- fwb:::make_gen_weights(wt)

    set.seed(57, "L'Ecuyer-CMRG")
    a <- gen(n, 10L, strata)

    set.seed(57, "L'Ecuyer-CMRG")
    b <- gen(n, 10L, NULL)

    expect_identical(a, b, info = wt)
  }
})

test_that("batch and per-replicate weight generation draw the same stream", {
  #`simple = TRUE` calls the generator `R` times with `R = 1`; `simple = FALSE`
  #calls it once with `R = R`. `fwb.array()` assumes these two consume the RNG
  #identically -- that is why the continuous generators fill their matrices with
  #`byrow = TRUE`. Losing that would corrupt `fwb.array()` and BCa intervals
  #without any visible error.
  n <- 25L
  R <- 8L

  for (wt in c("exp", "poisson", "mammen", "beta", "power")) {
    gen <- fwb:::make_gen_weights(wt)

    set.seed(58, "L'Ecuyer-CMRG")
    batch <- gen(n, R, NULL)

    set.seed(58, "L'Ecuyer-CMRG")
    seq_rows <- do.call("rbind", lapply(seq_len(R), function(i) drop(gen(n, 1L, NULL))))

    expect_equal(unname(batch), unname(seq_rows),
                 tolerance = fwb_eps(), info = wt)
  }
})

test_that("`\"multinom\"` batch and per-replicate generation deliberately differ", {
  #`"multinom"` assigns into its index matrix with `dim() <-` (column-major) rather than
  #`byrow = TRUE`, so the batch's first row is *not* the first `n` draws. That is not an
  #oversight: it is a verbatim transcription of `boot:::ordinary.array()`, and it is what
  #makes `wtype = "multinom"` reproduce `boot::boot(., stype = "f")` exactly
  #(`test-boot.R`). The cost is that batch and per-replicate generation disagree, which
  #`fwb.array()` handles by never batching for this `wtype`.
  n <- 25L
  R <- 8L

  gen <- fwb:::make_gen_weights("multinom")

  set.seed(59, "L'Ecuyer-CMRG")
  batch <- gen(n, R, NULL)

  set.seed(59, "L'Ecuyer-CMRG")
  seq_rows <- do.call("rbind", lapply(seq_len(R), function(i) drop(gen(n, 1L, NULL))))

  expect_not_equal(unname(batch), unname(seq_rows))

  #Both are valid multinomial weight matrices -- each row still sums to `n`. That is why
  #the mismatch is invisible in any marginal summary, and why it went unnoticed.
  expect_equal(rowSums(batch), rep(n, R))
  expect_equal(rowSums(seq_rows), rep(n, R))

  #And the batch order is `boot`'s, which is the property that must not change.
  skip_if_not_installed("boot")

  set.seed(59, "L'Ecuyer-CMRG")
  b_idx <- boot:::ordinary.array(n, R, rep(1, n))

  expect_identical(unname(batch), unname(t(apply(b_idx, 1L, tabulate, n))))
})

test_that("`tabulate_rows()` matches a row-by-row tabulation on both branches", {
  #`tabulate_rows()` picks between two counting strategies at `n = 256`. Both must
  #agree with the obvious implementation, so sizes here straddle the threshold.
  reference <- function(i, n) {
    out <- matrix(0L, nrow = nrow(i), ncol = n)
    for (r in seq_len(nrow(i))) {
      out[r, ] <- tabulate(i[r, ], n)
    }
    out
  }

  for (n in c(2L, 37L, 255L, 256L, 257L, 400L)) {
    for (R in c(1L, 3L, 40L)) {
      set.seed(n + R, "L'Ecuyer-CMRG")
      i <- sample.int(n, n * R, replace = TRUE)
      dim(i) <- c(R, n)

      expect_identical(fwb:::tabulate_rows(i, n), reference(i, n),
                       info = sprintf("n = %d, R = %d", n, R))
    }
  }
})

test_that("`\"multinom\"` weights are `R x n` even when `n` is 1", {
  #Regression test. The tabulation used to be `t(apply(i, 1L, tabulate, n))`, and
  #`apply()` simplifies a length-1 result to a vector, so at `n = 1` the transpose
  #produced a `1 x R` matrix instead of `R x 1`. Every other `wtype` returned
  #`R x 1`, and `fwb()` failed with "subscript out of bounds" on a one-row dataset
  #for `"multinom"` alone.
  for (wt in fwb:::.w_types()) {
    gen <- fwb:::make_gen_weights(wt)

    set.seed(71, "L'Ecuyer-CMRG")
    w <- gen(1L, 5L, NULL)

    expect_identical(dim(w), c(5L, 1L), info = wt)
  }

  #`fwb()` itself now refuses a one-unit dataset before any weight is drawn, so the
  #shape above is only reachable through the generator. The check stays because
  #`tabulate_rows()` is what makes it right, and nothing else guarantees it.
  expect_err(fwb(data.frame(x = 5, y = 2),
                 function(data, w) c(m = w_mean(data$x, w)),
                 R = 4L, wtype = "multinom", verbose = FALSE),
             "must have more than one unit")
})

test_that("`\"multinom\"` matches `boot` on both `tabulate_rows()` branches", {
  #`test-boot.R` and `test-fwb.R` check parity on a 2000-row fixture, which only ever
  #exercises the `n > 256` branch. Parity at small `n` needs its own test, or the
  #single-`tabulate()` branch is unverified against `boot`.
  skip_on_cran()
  skip_if_not_installed("boot")

  boot_fun <- function(data, w) {
    coef(lm(y ~ x1 + x2, data = data, weights = w))
  }

  #60 rows takes the first branch, 300 the second.
  for (n in c(60L, 300L)) {
    d <- fwb_test_df(n)

    set.seed(1234, "L")
    f <- fwb(d, boot_fun, R = 12L, verbose = FALSE, wtype = "multinom")

    set.seed(1234, "L")
    b <- boot::boot(d, boot_fun, R = 12L, stype = "f")

    expect_equal(f[["t"]], b[["t"]], tolerance = fwb_eps(),
                 ignore_attr = TRUE, info = sprintf("n = %d", n))

    #Same, stratified.
    set.seed(1234, "L")
    fs <- fwb(d, boot_fun, R = 12L, verbose = FALSE, wtype = "multinom",
              strata = s)

    set.seed(1234, "L")
    bs <- boot::boot(d, boot_fun, R = 12L, stype = "f", strata = d[["s"]])

    expect_equal(fs[["t"]], bs[["t"]], tolerance = fwb_eps(),
                 ignore_attr = TRUE, info = sprintf("n = %d, stratified", n))
  }
})

test_that("`get_fwb_wtype()` and `set_fwb_wtype()` round-trip", {
  op <- set_fwb_wtype("beta")
  defer_call(quote(options(.)), op)

  expect_identical(get_fwb_wtype(), "beta")
  expect_identical(getOption("fwb_wtype"), "beta")

  #Abbreviations are documented as allowed.
  invisible(set_fwb_wtype("mam"))
  expect_identical(get_fwb_wtype(), "mammen")

  d <- fwb_test_df()

  set.seed(60, "L'Ecuyer-CMRG")
  f <- fwb(d, lm_stat, R = 8, verbose = FALSE)

  expect_identical(f[["wtype"]], "mammen")
  expect_identical(get_fwb_wtype(f), "mammen")

  #An explicit `wtype` still wins over the option.
  set.seed(60, "L'Ecuyer-CMRG")
  f2 <- fwb(d, lm_stat, R = 8, verbose = FALSE, wtype = "exp")

  expect_identical(f2[["wtype"]], "exp")
  expect_identical(get_fwb_wtype(f2), "exp")

  #`set_fwb_wtype()` returns the previous value so `options()` restores it.
  options(op)
  expect_identical(get_fwb_wtype(), "exp")

  expect_err(set_fwb_wtype("nope"), "should be one of")
})

test_that("`vcovFWB()` honors the `fwb_wtype` option", {
  d <- fwb_test_df()
  fit <- lm(y ~ x1 + x2, data = d)

  op <- set_fwb_wtype("multinom")
  defer_call(quote(options(.)), op)

  set.seed(61, "L'Ecuyer-CMRG")
  v_opt <- vcovFWB(fit, R = 20)

  options(op)

  set.seed(61, "L'Ecuyer-CMRG")
  v_arg <- vcovFWB(fit, R = 20, wtype = "multinom")

  expect_equal(v_opt, v_arg, tolerance = fwb_eps())

  set.seed(61, "L'Ecuyer-CMRG")
  v_exp <- vcovFWB(fit, R = 20)

  expect_not_equal(v_opt, v_exp)
})

test_that("the multinomial bootstrap is enumerated when `R` reaches the number of samples", {
  #The multinomial weights have finitely many possible values, so at small `n` there is
  #no reason to sample: `fwb()` enumerates them and reduces `R` to that count.
  st <- function(data, w) c(m = w_mean(data[["x"]], w))

  d <- data.frame(x = c(1, 3))

  expect_message({
    set.seed(1, "L'Ecuyer-CMRG")
    f <- fwb(d, st, R = 999L, wtype = "multinom", verbose = FALSE)
  }, "only 4 distinct")

  expect_identical(f[["R"]], 4L)
  expect_identical(nrow(f[["t"]]), 4L)
  expect_true(isTRUE(attr(f, "exhaustive")))

  #`simple` is meaningless without randomness, and is reported as `FALSE`.
  expect_false(isTRUE(attr(f, "simple")))

  #The four samples of two units, in a fixed order. Two of them are `(1, 1)`: the
  #enumeration is over ordered resamples, which is what reproduces the multinomial
  #distribution. Keeping only the three distinct vectors would not.
  expect_identical(fwb.array(f),
                   matrix(c(2L, 1L, 1L, 0L,
                            0L, 1L, 1L, 2L), nrow = 4L))

  #No random numbers are drawn, so the result cannot depend on the seed.
  set.seed(987654, "L'Ecuyer-CMRG")
  f2 <- suppressMessages(fwb(d, st, R = 999L, wtype = "multinom", verbose = FALSE))

  expect_identical(f[["t"]], f2[["t"]])
  expect_identical(fwb.array(f), fwb.array(f2))

  #Not triggered when `R` is below the count, and never for any other weight type.
  set.seed(2, "L'Ecuyer-CMRG")
  f3 <- fwb(data.frame(x = 1:3), st, R = 10L, wtype = "multinom", verbose = FALSE)

  expect_identical(f3[["R"]], 10L)
  expect_false(isTRUE(attr(f3, "exhaustive")))

  #A statistic that survives an all-zero weight vector, which `"poisson"` can produce
  #at this size; `w_mean()` would return `NA` and warn about it.
  sum_st <- function(data, w) c(s = sum(w * data[["x"]]))

  for (wt in setdiff(fwb:::.w_types(), "multinom")) {
    set.seed(3, "L'Ecuyer-CMRG")
    fw <- fwb(data.frame(x = 1:3), sum_st, R = 40L, wtype = wt, verbose = FALSE)

    expect_identical(fw[["R"]], 40L, info = wt)
    expect_false(isTRUE(attr(fw, "exhaustive")), info = wt)
  }
})

test_that("the enumeration is complete and correctly weighted", {
  st <- function(data, w) c(m = w_mean(data[["x"]], w))

  #With `n = 3` there are 27 ordered resamples and 10 distinct weight vectors. Each
  #vector must appear as often as its multinomial coefficient: 1 for `(3,0,0)`-type,
  #3 for `(2,1,0)`-type, 6 for `(1,1,1)`.
  d <- data.frame(x = c(1, 2, 4))

  f <- suppressMessages(fwb(d, st, R = 999L, wtype = "multinom", verbose = FALSE))

  a <- fwb.array(f)

  expect_identical(f[["R"]], 27L)
  expect_identical(nrow(a), 27L)
  expect_true(all(rowSums(a) == 3L))

  counts <- table(apply(a, 1L, paste, collapse = ""))

  expect_identical(sum(counts), 27L)
  expect_length(counts, 10L)
  expect_identical(unname(counts[["111"]]), 6L)
  expect_identical(unname(counts[["300"]]), 1L)
  expect_identical(unname(counts[["210"]]), 3L)

  #Strata multiply: the enumeration is the product of the per-stratum enumerations,
  #and each stratum's weights still sum to its size in every replicate.
  ds <- data.frame(x = c(1, 2, 3, 4), s = factor(c("a", "a", "b", "b")))

  fs <- suppressMessages(fwb(ds, st, R = 999L, wtype = "multinom", strata = s,
                             verbose = FALSE))

  as <- fwb.array(fs)

  expect_identical(fs[["R"]], 16L)
  expect_true(all(rowSums(as[, 1:2]) == 2L))
  expect_true(all(rowSums(as[, 3:4]) == 2L))

  #Clusters enumerate at the cluster scale, and members of a cluster share a weight.
  dc <- data.frame(x = 1:6, g = factor(rep(1:3, each = 2L)))

  fc <- suppressMessages(fwb(dc, st, R = 999L, wtype = "multinom", cluster = g,
                             verbose = FALSE))

  ac <- fwb.array(fc)

  expect_identical(fc[["R"]], 27L)
  expect_true(all(ac[, 1] == ac[, 2] & ac[, 3] == ac[, 4] & ac[, 5] == ac[, 6]))
})

test_that("`n_multinom_resamples()` counts what it claims to", {
  for (n in 2:7) {
    expect_identical(fwb:::n_multinom_resamples(n), as.double(n)^n)
  }

  #Independent across strata, so the count is the product.
  s <- factor(c("a", "a", "b", "b", "b"))

  expect_identical(fwb:::n_multinom_resamples(5L, s), 2^2 * 3^3)

  #A single level is the same as no strata at all.
  expect_identical(fwb:::n_multinom_resamples(4L, factor(rep("a", 4L))),
                   fwb:::n_multinom_resamples(4L))

  #Enumerating actually produces that many rows.
  expect_identical(nrow(fwb:::gen_all_multinom_weights(4L)), 256L)
  expect_identical(nrow(fwb:::gen_all_multinom_weights(5L, s)), 108L)
})

test_that("every stratum of one leaves nothing to bootstrap", {
  #Each stratum's only possible resample is itself, so all weights are 1 in every
  #replicate. Enumerating would give a single replicate; refusing is more honest.
  d <- data.frame(x = 1:3, s = factor(1:3))

  expect_err(fwb(d, function(data, w) c(m = w_mean(data[["x"]], w)),
                 R = 99L, wtype = "multinom", strata = s, verbose = FALSE),
             "only one possible multinomial bootstrap sample")
})
