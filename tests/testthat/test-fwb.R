test_that("fwb() works", {
  eps <- if (capabilities("long.double")) 1e-8 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  test_data$clus <- sample(1:50, nrow(test_data), replace = TRUE)

  boot_fun <- function(data, w = NULL) {
    fit <- glm(Y_B ~ A + X1 + X2 + X3 + X4, data = data,
               family = quasibinomial("probit"), weights = w)
    coef(fit)
  }

  set.seed(1234, "L")
  expect_no_condition({
    f0 <- fwb(test_data, boot_fun, R = 100, verbose = FALSE)
  })

  expect_identical(names(f0),
                   c("t0", "t", "R", "data", "statistic", "call", "cluster",
                     "strata", "wtype"))

  expect_equal(length(f0[["t0"]]), length(boot_fun(test_data)))

  expect_equal(ncol(f0[["t"]]), length(f0[["t0"]]))
  expect_equal(nrow(f0[["t"]]), f0[["R"]])
  expect_equal(f0[["data"]], test_data)
  expect_equal(f0[["statistic"]], boot_fun)
  expect_null(f0[["cluster"]])
  expect_null(f0[["strata"]])
  expect_equal(f0[["wtype"]], "exp")
  expect_true(attr(f0, "simple", TRUE))

  set.seed(1234, "L")
  expect_no_condition({
    f1 <- fwb(test_data, boot_fun, R = 100, verbose = FALSE,
              simple = FALSE)
  })

  expect_false(attr(f1, "simple", TRUE))

  #`simple` now changes the replicates. `simple = FALSE` draws the whole weight matrix
  #in one call (which is what keeps `wtype = "multinom"` identical to `boot`);
  #`simple = TRUE` draws each replicate from its own recorded stream, so that the
  #results do not depend on the backend. The two consume the stream differently, so
  #they no longer coincide -- but each is reproducible on its own, and both are valid
  #draws from the same weight distribution.
  expect_not_equal(f1[["t"]], f0[["t"]], tolerance = eps)

  expect_equal(f1[["t0"]], f0[["t0"]], tolerance = eps)
  expect_equal(dim(f1[["t"]]), dim(f0[["t"]]))
  expect_equal(colMeans(f1[["t"]]), colMeans(f0[["t"]]), tolerance = 1e-1)

  set.seed(1234, "L")
  expect_no_condition({
    f2 <- fwb(test_data, function(data, w) c(boot_fun(data, w), w),
              R = 100, verbose = FALSE, cluster = clus)
  })

  expect_identical(names(f2),
                   c("t0", "t", "R", "data", "statistic", "call", "cluster",
                     "strata", "wtype"))

  expect_equal(length(f2[["t0"]]), length(boot_fun(test_data)) + nrow(test_data))

  expect_equal(ncol(f2[["t"]]), length(f2[["t0"]]))
  expect_equal(nrow(f2[["t"]]), f2[["R"]])
  expect_equal(f2[["data"]], test_data)
  expect_failure(expect_null(f2[["cluster"]]))
  expect_null(f2[["strata"]])
  expect_equal(f2[["wtype"]], "exp")
  expect_true(attr(f2, "simple", TRUE))

  #Test that weights in each cluster are the same
  expect_true(all(apply(f2$t[,-(1:6)], 1, tapply, f2$cluster, function(z) all(z == z[1L]))))

  set.seed(1234, "L")
  expect_no_condition({
    f3 <- fwb(test_data, boot_fun, R = 100, verbose = FALSE,
              wtype = "mammen")
  })

  expect_not_equal(f0$t, f3$t, tolerance = eps,
                   ignore_attr = TRUE)
})

test_that("parallel works", {
  skip_on_cran()
  eps <- if (capabilities("long.double")) 1e-8 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  boot_fun <- function(data, w = NULL) {
    fit <- glm(Y_B ~ A + X1 + X2 + X3 + X4, data = data,
               family = quasibinomial("probit"), weights = w)
    coef(fit)
  }

  #Weights are drawn from streams recorded before any work is dispatched, so `cl` has
  #no bearing on the result: every backend and worker count gives the same replicates
  #as sequential evaluation, for either value of `simple`. `test-replicability.R`
  #covers the full grid; this is the smoke test.
  for (simple in c(TRUE, FALSE)) {
    set.seed(1234, "L")
    ref <- fwb(test_data, boot_fun, R = 100, verbose = FALSE, simple = simple)

    set.seed(1234, "L")
    expect_no_condition({
      f_int <- fwb(test_data, boot_fun, R = 100, verbose = FALSE,
                   cl = 2, simple = simple)
    })

    expect_equal(f_int[["t"]], ref[["t"]], tolerance = eps)

    cl <- parallel::makeCluster(2)

    set.seed(1234, "L")
    expect_no_condition({
      f_clus <- fwb(test_data, boot_fun, R = 100, verbose = FALSE,
                    cl = cl, simple = simple)
    })

    parallel::stopCluster(cl)

    expect_equal(f_clus[["t"]], ref[["t"]], tolerance = eps)
  }

  #`set.seed()` in the calling session is what makes a `cluster` run reproducible;
  #`parallel::clusterSetRNGStream()` no longer has any effect, because the workers'
  #own streams are never used to draw weights.
  cl <- parallel::makeCluster(2)
  on.exit(parallel::stopCluster(cl))

  parallel::clusterSetRNGStream(cl, 1234)
  a <- fwb(test_data, boot_fun, R = 100, verbose = FALSE, cl = cl)

  parallel::clusterSetRNGStream(cl, 1234)
  b <- fwb(test_data, boot_fun, R = 100, verbose = FALSE, cl = cl)

  expect_not_equal(a[["t"]], b[["t"]], tolerance = eps)

  set.seed(99, "L")
  c1 <- fwb(test_data, boot_fun, R = 100, verbose = FALSE, cl = cl)

  set.seed(99, "L")
  c2 <- fwb(test_data, boot_fun, R = 100, verbose = FALSE, cl = cl)

  expect_equal(c1[["t"]], c2[["t"]], tolerance = eps)
})

test_that("wtype = 'multinom' replcates boot::boot()", {
  skip_on_cran()
  eps <- if (capabilities("long.double")) 1e-8 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(123, "L")

  clus <- sample(1:50, nrow(test_data), replace = TRUE)

  boot_fun <- function(data, w) {
    fit <- glm(Y_B ~ A + X1 + X2 + X3 + X4, data = data,
               family = quasibinomial("probit"), weights = w)
    coef(fit)
  }

  cl <- parallel::makeCluster(2)
  on.exit(parallel::stopCluster(cl))

  set.seed(1234, "L")
  expect_no_condition({
    f0 <- fwb(test_data, boot_fun, R = 10, verbose = FALSE,
             wtype = "multinom", simple = TRUE)
  })

  #Without strata
  set.seed(1234, "L")
  expect_no_condition({
    f <- fwb(test_data, boot_fun, R = 10, verbose = FALSE,
             wtype = "multinom")
  })

  expect_not_equal(f$t, f0$t)

  set.seed(1234, "L")
  expect_no_condition({
    b <- boot::boot(test_data, boot_fun, R = 10,
                    stype = "f")
  })

  expect_equal(f$t, b$t, tolerance = eps,
               ignore_attr = TRUE)

  #With strata
  set.seed(1234, "L")
  expect_no_condition({
    f <- fwb(test_data, boot_fun, R = 10, verbose = FALSE,
             wtype = "multinom", strata = A)
  })

  set.seed(1234, "L")
  expect_no_condition({
    b <- boot::boot(test_data, boot_fun, R = 10,
                    stype = "f", strata = test_data$A)
  })

  expect_equal(f$t, b$t, tolerance = eps,
               ignore_attr = TRUE)
})

test_that("drop0 works as expected", {
  skip_on_cran()
  eps <- if (capabilities("long.double")) 1e-8 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(123, "L")

  clus <- sample(1:50, nrow(test_data), replace = TRUE)

  boot_fun <- function(data, w, abort = TRUE) {
    if (abort && any(w[!is.na(w)] == 0)) {
      stop("bad w")
    }

    fit <- glm(Y_B ~ A + X1 + X2 + X3 + X4, data = data,
               family = quasibinomial("probit"), weights = w)
    coef(fit)
  }

  set.seed(1234, "L")
  expect_error({
    f <- fwb(test_data, boot_fun, R = 20, verbose = FALSE,
             wtype = "multinom")
  }, .w("bad w"))

  set.seed(1234, "L")
  expect_no_condition({
    fF <- fwb(test_data, boot_fun, R = 20, verbose = FALSE,
              wtype = "multinom", drop0 = FALSE, abort = FALSE)
  })

  set.seed(1234, "L")
  expect_no_condition({
    fT <- fwb(test_data, boot_fun, R = 20, verbose = FALSE,
             wtype = "multinom", drop0 = TRUE)
  })

  expect_equal(fT$t, fF$t, tolerance = eps,
               ignore_attr = TRUE)

  set.seed(1234, "L")
  expect_no_condition({
    fNA <- fwb(test_data, boot_fun, R = 20, verbose = FALSE,
              wtype = "multinom", drop0 = NA)
  })

  expect_equal(fT$t, fNA$t, tolerance = eps,
               ignore_attr = TRUE)
})
