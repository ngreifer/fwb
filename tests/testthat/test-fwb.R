test_that("fwb() works", {
  eps <- if (capabilities("long.double")) 1e-8 else 1e-1

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
                   c("t0", "t", "R", "data", "seed", "statistic", "call", "cluster",
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

  expect_equal(f1[-7], f0[-7], tolerance = eps)
  expect_false(attr(f1, "simple", TRUE))

  set.seed(1234, "L")
  expect_no_condition({
    f2 <- fwb(test_data, function(data, w) c(boot_fun(data, w), w),
              R = 100, verbose = FALSE, cluster = clus)
  })

  expect_identical(names(f2),
                   c("t0", "t", "R", "data", "seed", "statistic", "call", "cluster",
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
  eps <- if (capabilities("long.double")) 1e-8 else 1e-1

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

  set.seed(1234, "L")
  expect_no_condition({
    f1 <- fwb(test_data, boot_fun, R = 100, verbose = FALSE,
              cl = 2, simple = FALSE)
  })

  expect_equal(f1[-7], f0[-7], tolerance = eps)
  expect_false(attr(f1, "simple", TRUE))

  #Using cl = int
  set.seed(1234, "L")
  expect_no_condition({
    f2 <- fwb(test_data, boot_fun, R = 100, verbose = FALSE,
              cl = 2, simple = TRUE)
  })

  set.seed(1234, "L")
  expect_no_condition({
    f3 <- fwb(test_data, boot_fun, R = 100, verbose = FALSE,
              cl = 2, simple = TRUE)
  })

  expect_equal(f2, f3, tolerance = eps)

  #Using a cluster
  cl <- parallel::makeCluster(2)
  on.exit(parallel::stopCluster(cl))

  parallel::clusterSetRNGStream(cl, 1234)
  expect_no_condition({
    f2 <- fwb(test_data, boot_fun, R = 100, verbose = FALSE,
              cl = cl, simple = TRUE)
  })

  parallel::clusterSetRNGStream(cl, 1234)
  expect_no_condition({
    f3 <- fwb(test_data, boot_fun, R = 100, verbose = FALSE,
              cl = cl, simple = TRUE)
  })

  expect_equal(f2, f3, tolerance = eps)
})

test_that("wtype = 'multinom' replcates boot::boot()", {
  skip_on_cran()
  eps <- if (capabilities("long.double")) 1e-8 else 1e-1

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
  eps <- if (capabilities("long.double")) 1e-8 else 1e-1

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
