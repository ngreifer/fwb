test_that("`w_*()` functions work like base R weighted functions", {
  set.seed(123456)
  eps <- if (capabilities("long.double")) 1e-8 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  sw <- test_data$SW
  w1 <- rep(1, length(sw))
  x <- test_data$X9

  # w_mean()
  expect_equal(w_mean(x), mean(x),
               tolerance = eps)

  expect_equal(w_mean(x, sw), weighted.mean(x, sw),
               tolerance = eps)

  expect_equal(w_mean(x), w_mean(x, w1),
               tolerance = eps)

  # w_var()
  expect_equal(w_var(x), var(x),
               tolerance = eps)

  expect_equal(w_var(x, sw), cov.wt(as.matrix(x), sw)$cov[1],
               tolerance = eps)

  expect_equal(w_var(x), w_var(x, w1),
               tolerance = eps)

  # w_sd()
  expect_equal(w_sd(x), sd(x),
               tolerance = eps)

  expect_equal(w_sd(x, sw), sqrt(w_var(x, sw)),
               tolerance = eps)

  expect_equal(w_sd(x), w_sd(x, w1),
               tolerance = eps)

  # w_cov()
  X <- test_data[4:8]

  expect_equal(w_cov(X),
               cov(X),
               tolerance = eps)

  expect_equal(w_cov(X, sw),
               cov.wt(X, sw)$cov,
               tolerance = eps)

  expect_equal(w_cov(X),
               w_cov(X, w1),
               tolerance = eps)

  # w_cor()
  expect_equal(w_cor(X),
               cor(X),
               tolerance = eps)

  expect_equal(w_cor(X, sw),
               cov.wt(X, sw, cor = TRUE)$cor,
               tolerance = eps)

  expect_equal(w_cor(X, sw), cov2cor(w_cov(X, sw)),
               tolerance = eps)

  expect_equal(w_cor(X), w_cor(X, w1),
               tolerance = eps)

  # w_quantile()
  expect_equal(w_quantile(x), quantile(x),
               tolerance = eps)

  expect_equal(w_quantile(x), w_quantile(x, w1),
               tolerance = eps)

  if (rlang::is_installed("ggdist")) {
    expect_equal(w_quantile(x, sw),
                 ggdist::weighted_quantile(x, weights = sw),
                 tolerance = eps)
  }

  # w_median()
  expect_equal(w_median(x), median(x),
               tolerance = eps)

  expect_equal(w_median(x), w_median(x, w1),
               tolerance = eps)

  expect_equal(w_median(x, sw),
               w_quantile(x, sw, .5, names = FALSE),
               tolerance = eps)

  # w_std()
  expect_equal(w_std(x),
               drop(scale(x)),
               tolerance = eps,
               ignore_attr = TRUE)

  expect_equal(w_std(x),
               w_std(x, w1),
               tolerance = eps)

  expect_equal(w_std(x, sw),
               drop(scale(x, center = weighted.mean(x, sw),
                          scale = sqrt(cov.wt(as.matrix(x), sw)$cov[1]))),
               tolerance = eps,
               ignore_attr = TRUE)

  expect_equal(w_std(x, sw),
               (x - w_mean(x, sw)) / w_sd(x, sw),
               tolerance = eps)
})

test_that("`w_*()` functions automatically capture weights inside `fwb()`", {
  set.seed(12345)
  eps <- if (capabilities("long.double")) 1e-8 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  sw <- test_data$SW
  w1 <- rep(1, length(sw))

  # w_mean()
  bootfun <- function(data, .w = rep(1, nrow(data))) {
    c(mean = mean(data$X9),
      weighted.mean = weighted.mean(data$X9, .w),
      w_mean_w = w_mean(data$X9, .w),
      w_mean_1 = w_mean(data$X9, w1),
      w_mean_auto = w_mean(data$X9))
  }

  expect_no_condition({
    boot_out <- fwb(test_data, bootfun, R = 20, verbose = FALSE)
  })

  expect_equal(boot_out$t0,
               rep(mean(test_data$X9), length(boot_out$t0)),
               tolerance = eps,
               ignore_attr = TRUE)

  expect_equal(boot_out$t[,"mean"],
               boot_out$t[,"w_mean_1"],
               tolerance = eps)

  expect_equal(boot_out$t[,"weighted.mean"],
               boot_out$t[,"w_mean_w"],
               tolerance = eps)

  expect_equal(boot_out$t[,"weighted.mean"],
               boot_out$t[,"w_mean_auto"],
               tolerance = eps)

  expect_true(sd(boot_out$t[,"w_mean_auto"]) > 0)

  expect_equal(sd(boot_out$t[,"mean"]),
               0,
               tolerance = eps)

  # w_var()
  bootfun <- function(data, .w = rep(1, nrow(data))) {
    c(var = var(data$X9),
      cov.wt = cov.wt(as.matrix(data$X9), .w)$cov[1],
      w_var_w = w_var(data$X9, .w),
      w_var_1 = w_var(data$X9, w1),
      w_var_auto = w_var(data$X9))
  }

  expect_no_condition({
    boot_out <- fwb(test_data, bootfun, R = 20, verbose = FALSE)
  })

  expect_equal(boot_out$t0,
               rep(var(test_data$X9), length(boot_out$t0)),
               tolerance = eps,
               ignore_attr = TRUE)

  expect_equal(boot_out$t[,"var"],
               boot_out$t[,"w_var_1"],
               tolerance = eps)

  expect_equal(boot_out$t[,"cov.wt"],
               boot_out$t[,"w_var_w"],
               tolerance = eps)

  expect_equal(boot_out$t[,"w_var_w"],
               boot_out$t[,"w_var_auto"],
               tolerance = eps)

  expect_true(sd(boot_out$t[,"w_var_auto"]) > 0)

  expect_equal(sd(boot_out$t[,"var"]),
               0,
               tolerance = eps)

  # w_median()
  bootfun <- function(data, .w = rep(1, nrow(data))) {
    c(median = median(data$X9),
      w_median_1 = w_median(data$X9, w1),
      w_median_w = w_median(data$X9, .w),
      w_median_auto = w_median(data$X9))
  }

  expect_no_condition({
    boot_out <- fwb(test_data, bootfun, R = 20, verbose = FALSE)
  })

  expect_equal(boot_out$t0,
               rep(median(test_data$X9), length(boot_out$t0)),
               tolerance = eps,
               ignore_attr = TRUE)

  expect_equal(boot_out$t[,"median"],
               boot_out$t[,"w_median_1"],
               tolerance = eps)

  expect_equal(boot_out$t[,"w_median_w"],
               boot_out$t[,"w_median_auto"],
               tolerance = eps)

  expect_true(sd(boot_out$t[,"w_median_auto"]) > 0)

  expect_equal(sd(boot_out$t[,"median"]),
               0,
               tolerance = eps)

  # w_std()
  bootfun <- function(data, .w = rep(1, nrow(data))) {
    fit1 <- lm(Y_C ~ A * w_std(X1), data = test_data, weights = .w)
    fit2 <- lm(Y_C ~ A * w_std(X1, .w), data = test_data, weights = .w)
    fit3 <- lm(Y_C ~ A * w_std(X1, w1), data = test_data, weights = .w)
    fit4 <- lm(Y_C ~ A * X1_scale, data = transform(test_data,
                                                    X1_scale = drop(scale(X1, center = weighted.mean(X1, .w),
                                                                          scale = sqrt(cov.wt(as.matrix(X1), .w)$cov[1])))),
               weights = .w)
    fit5 <- lm(Y_C ~ A * scale(X1), data = test_data, weights = .w)

    c(auto = unname(coef(fit1)["A"]),
      manual = unname(coef(fit2)["A"]),
      uniform = unname(coef(fit3)["A"]),
      transform = unname(coef(fit4)["A"]),
      ignore = unname(coef(fit5)["A"]))
  }

  set.seed(111)
  expect_no_condition({
    boot_out <- fwb(test_data, bootfun, R = 20, verbose = FALSE)
  })

  expect_equal(boot_out$t0,
               rep(coef(lm(Y_C ~ A * scale(X1), data = test_data))["A"], length(boot_out$t0)),
               tolerance = eps,
               ignore_attr = TRUE)

  expect_equal(boot_out$t[,"auto"],
               boot_out$t[,"manual"],
               tolerance = eps)

  expect_equal(boot_out$t[,"auto"],
               boot_out$t[,"transform"],
               tolerance = eps)

  expect_equal(boot_out$t[,"uniform"],
               boot_out$t[,"ignore"],
               tolerance = eps)

  expect_not_equal(boot_out$t[,"auto"],
                   boot_out$t[,"uniform"])

  fit <- lm(Y_C ~ A * w_std(X1), data = test_data)
  set.seed(111)
  v <- vcovFWB(fit, R = 20)

  expect_equal(sqrt(v["A","A"]),
               sd(boot_out$t[, "auto"]),
               tolerance = eps)

  expect_not_equal(sqrt(v["A","A"]),
                   sd(boot_out$t[, "uniform"]))
})
