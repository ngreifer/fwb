#`plot.fwb()`.
#
#The plot cannot be checked for what it looks like without a snapshot, and a base-
#graphics snapshot is not worth its maintenance cost here. What is worth pinning is
#that every combination of arguments draws without erroring, restores `par()`, and
#refuses the inputs it cannot draw -- the failure modes that would otherwise show up
#as a broken vignette or a leftover `Rplots.pdf`.

test_that("`plot.fwb()` draws both panels and restores `par()`", {
  local_null_device()

  d <- fwb_test_df()

  set.seed(103, "L'Ecuyer-CMRG")
  f <- fwb(d, lm_stat, R = 200L, verbose = FALSE, simple = FALSE)

  before <- graphics::par("mfrow")

  expect_no_error({
    res <- plot(f, index = 2L)
  })

  #`plot()` sets `mfrow` for the side-by-side layout and must put it back.
  expect_identical(graphics::par("mfrow"), before)

  #The `<fwb>` object is returned invisibly.
  expect_identical(res, f)
})

test_that("`plot.fwb()` honors `type`", {
  local_null_device()

  d <- fwb_test_df()

  set.seed(103, "L'Ecuyer-CMRG")
  f <- fwb(d, lm_stat, R = 200L, verbose = FALSE, simple = FALSE)

  for (ty in list("hist", "qq", c("hist", "qq"))) {
    expect_no_error({
      plot(f, index = 2L, type = ty)
    })
  }

  expect_err(plot(f, index = 2L, type = "nope"), "should be")
})

test_that("`plot.fwb()` honors `qdist` and `df`", {
  local_null_device()

  d <- fwb_test_df()

  #A positive statistic so that a chi-squared Q-Q plot is meaningful.
  chisq_stat <- function(data, w) {
    c(v = w_var(data[["y"]], w))
  }

  set.seed(104, "L'Ecuyer-CMRG")
  f <- fwb(d, chisq_stat, R = 200L, verbose = FALSE, simple = FALSE)

  expect_no_error({
    plot(f, type = "qq", qdist = "norm")
  })

  #With `df` supplied, and with `df` estimated by maximum likelihood.
  expect_no_error({
    plot(f, type = "qq", qdist = "chisq", df = 3)
  })

  expect_no_error({
    plot(f, type = "qq", qdist = "chisq")
  })

  expect_err(plot(f, type = "qq", qdist = "nope"), "should be one of")
  expect_err(plot(f, type = "qq", qdist = "chisq", df = "three"),
             "must be a number")
})

test_that("`plot.fwb()` honors `nclass`", {
  local_null_device()

  d <- fwb_test_df()

  set.seed(103, "L'Ecuyer-CMRG")
  f <- fwb(d, lm_stat, R = 200L, verbose = FALSE, simple = FALSE)

  for (nc in c(10L, 30L, 100L)) {
    expect_no_error({
      plot(f, index = 2L, type = "hist", nclass = nc)
    })
  }

  expect_err(plot(f, index = 2L, type = "hist", nclass = -5L), "must be")
  expect_err(plot(f, index = 2L, type = "hist", nclass = 2.5), "must be")
})

test_that("`plot.fwb()` selects one statistic at a time", {
  local_null_device()

  d <- fwb_test_df()

  set.seed(103, "L'Ecuyer-CMRG")
  f <- fwb(d, lm_stat, R = 200L, verbose = FALSE, simple = FALSE)

  expect_no_error({
    plot(f, index = "x1")
  })

  expect_err(plot(f, index = 1:2), "`index` must have length one")
  expect_err(plot(f, index = 99L), "`index` must be between 1 and")
  expect_err(plot(f, index = "nope"), "all entries in `index` must be the names")
})

test_that("`plot.fwb()` warns rather than drawing a degenerate distribution", {
  local_null_device()

  d <- fwb_test_df()

  set.seed(105, "L'Ecuyer-CMRG")
  f <- fwb(d, function(data, w) c(a = 1), R = 50L, verbose = FALSE)

  expect_wrn(res <- plot(f), "all values of")

  expect_identical(res, f)
})

test_that("`plot.fwb()` handles non-finite replicates", {
  local_null_device()

  d <- fwb_test_df()

  #Non-finite replicates are dropped before plotting, so the Q-Q panel needs one
  #theoretical quantile per *finite* replicate. It used to build them from `x$R`, leaving
  #the two vectors different lengths; `qqplot()` tolerates that by interpolating, so the
  #plot was drawn but its points were not the quantiles its axis implied.
  sometimes_inf <- function(data, w) {
    est <- coef(lm(y ~ x1, data = data, weights = w))[[2L]]

    c(v = if (w[[1L]] > 2.5) Inf else est)
  }

  set.seed(106, "L'Ecuyer-CMRG")
  f <- suppressWarnings(fwb(d, sometimes_inf, R = 200L, verbose = FALSE))

  n_finite <- sum(is.finite(f[["t"]]))

  expect_lt(n_finite, f[["R"]])

  #Intercept `qqplot()` to see what it is actually handed. Comparing the two lengths is
  #the whole assertion: equal to each other, and equal to the number of points there
  #genuinely are.
  seen <- NULL

  local_mocked_bindings(
    qqplot = function(x, y, ...) {
      seen <<- c(length(x), length(y))
      invisible(list(x = x, y = y))
    },
    qqline = function(...) invisible(NULL),
    .package = "fwb")

  expect_no_error({
    plot(f, type = c("hist", "qq"))
  })

  expect_identical(seen, c(n_finite, n_finite))
})
