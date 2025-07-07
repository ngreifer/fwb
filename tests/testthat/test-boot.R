test_that("fwb aligns with boot", {
  skip_if_not_installed("boot")

  eps <- if (capabilities("long.double")) 1e-8 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  set.seed(12345, "L")
  test_data <- test_data[sample.int(nrow(test_data), 500),]

  set_fwb_wtype("multinom")
  on.exit(set_fwb_wtype("exp"))

  boot_fun <- function(data, w = NULL) {
    fit <- glm(Y_B ~ A + X1 + X2 + X3 + X4, data = data,
               family = quasibinomial("logit"), weights = w)
    coef(fit)
  }

  set.seed(1234, "L")
  expect_no_condition({
    f0 <- fwb(test_data, boot_fun, R = 550, verbose = FALSE)
  })

  set.seed(1234, "L")
  expect_no_condition({
    b0 <- boot::boot(test_data, boot_fun, R = 550, stype = "f")
  })

  expect_equal(f0[["t"]], b0[["t"]],
               tolerance = eps, ignore_attr = TRUE)

  expect_equal(f0[["t0"]], b0[["t0"]],
               tolerance = eps, ignore_attr = TRUE)

  ci.types <- c("norm", "basic", "perc", "bca")

  expect_no_condition({
    f0.ci <- fwb.ci(f0, conf = .86, type = ci.types, index = 3)
  })

  expect_no_condition({
    b0.ci <- boot::boot.ci(b0, conf = .86, type = ci.types, index = 3)
  })

  expect_equal(f0.ci[-3], b0.ci[-3],
               tolerance = eps, ignore_attr = TRUE)

  expect_equal(get_ci(f0.ci), get_ci(b0.ci),
               tolerance = eps)

  # Test with hinv
  expect_no_condition({
    f0.ci <- fwb.ci(f0, conf = .86, type = ci.types, index = 3,
                    hinv = exp)
  })

  expect_no_condition({
    b0.ci <- boot::boot.ci(b0, conf = .86, type = ci.types, index = 3,
                           hinv = exp)
  })

  expect_equal(f0.ci[-3], b0.ci[-3],
               tolerance = eps, ignore_attr = TRUE)

  expect_equal(get_ci(f0.ci), get_ci(b0.ci),
               tolerance = eps)


  # Extreme endpoints
  for (i in c("basic", "perc", "bca")) {
    expect_warning({
      f0.ci <- fwb.ci(f0, conf = .999, type = i, index = 4)
    }, .w("xtreme order statistics used as endpoints"))

    expect_warning({
      b0.ci <- boot::boot.ci(b0, conf = .999, type = i, index = 4)
    }, .w("xtreme order statistics used as endpoints"))

    expect_equal(f0.ci[-3], b0.ci[-3],
                 tolerance = eps, ignore_attr = TRUE)

    expect_equal(get_ci(f0.ci), get_ci(b0.ci),
                 tolerance = eps, ignore_attr = TRUE)
  }

})
