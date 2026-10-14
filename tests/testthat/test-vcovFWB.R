test_that("vcovFWB() works", {
  eps <- fwb_eps()

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  test_data$clus <- sample(1:50, nrow(test_data), replace = TRUE)

  #vcovFWB produces same estimates as fwb(). Each replicate draws its weights from its
  #own stream, as in fwb(), so the draws do not reproduce those of sandwich::vcovBS().

  #lm()
  fit <- lm(Y_C ~ A + X1 + X2 + X3 + X4, data = test_data)

  boot_fun <- function(data, w = NULL) {
    boot_fit <- lm(Y_C ~ A + X1 + X2 + X3 + X4, data = data,
                   weights = w)
    coef(boot_fit)
  }

  # With wtype = "exp"
  set.seed(1234, "L")
  expect_no_condition({
    f0 <- fwb(test_data, boot_fun, R = 100, verbose = FALSE)
  })

  set.seed(1234, "L")
  expect_no_condition({
    v0 <- vcovFWB(fit, R = 100)
  })

  expect_equal(v0, vcov(f0), tolerance = eps)

  # With wtype = "multinom"
  set.seed(1234, "L")
  expect_no_condition({
    f0 <- fwb(test_data, boot_fun, R = 100, verbose = FALSE,
              wtype = "multinom", simple = TRUE)
  })

  set.seed(1234, "L")
  expect_no_condition({
    v0 <- vcovFWB(fit, R = 100,
                  wtype = "multinom")
  })

  expect_equal(v0, vcov(f0), tolerance = eps)


  # glm()
  fit <- glm(Y_B ~ A + X1 + X2 + X3 + X4, data = test_data,
             family = quasibinomial("probit"))

  boot_fun <- function(data, w = NULL) {
    boot_fit <- glm(Y_B ~ A + X1 + X2 + X3 + X4, data = data,
                    family = quasibinomial("probit"), weights = w)
    coef(boot_fit)
  }

  # With wtype = "exp"
  set.seed(1234, "L")
  expect_no_condition({
    f0 <- fwb(test_data, boot_fun, R = 100, verbose = FALSE)
  })

  set.seed(1234, "L")
  expect_no_condition({
    v0 <- vcovFWB(fit, R = 100)
  })

  expect_equal(v0, vcov(f0), tolerance = eps)

  # With wtype = "multinom"
  set.seed(1234, "L")
  expect_no_condition({
    f0 <- fwb(test_data, boot_fun, R = 100, verbose = FALSE,
              wtype = "multinom", simple = TRUE)
  })

  set.seed(1234, "L")
  expect_no_condition({
    v0 <- vcovFWB(fit, R = 100,
                  wtype = "multinom")
  })

  expect_equal(v0, vcov(f0), tolerance = eps)

  # coxph()
  library(survival)
  fit <- coxph(Surv(Y_S, Y_B) ~ A + X1 + X2 + X3 + X4, data = test_data)

  boot_fun <- function(data, w = NULL) {
    boot_fit <- coxph(Surv(Y_S, Y_B) ~ A + X1 + X2 + X3 + X4, data = data,
                      weights = w)
    coef(boot_fit)
  }

  # With wtype = "exp"
  set.seed(1234, "L")
  expect_no_condition({
    f0 <- fwb(test_data, boot_fun, R = 100, verbose = FALSE)
  })

  set.seed(1234, "L")
  expect_no_condition({
    v0 <- vcovFWB(fit, R = 100)
  })

  expect_equal(v0, vcov(f0), tolerance = eps)

  # set.seed(1234, "L")
  # expect_no_condition({
  #   vs <- sandwich::vcovBS(fit, R = 100, type = "fractional")
  # })
  #
  # expect_equal(v0, vs, tolerance = eps)

  # With wtype = "multinom"
  set.seed(1234, "L")
  expect_no_condition({
    f0 <- fwb(test_data, boot_fun, R = 100, verbose = FALSE,
              wtype = "multinom", simple = TRUE, drop0 = TRUE)
  })

  set.seed(1234, "L")
  expect_no_condition({
    v0 <- vcovFWB(fit, R = 100, drop0 = TRUE,
                  wtype = "multinom")
  })

  expect_equal(v0, vcov(f0), tolerance = eps)

  set.seed(1234, "L")
  expect_no_condition({
    fNA <- fwb(test_data, boot_fun, R = 100, verbose = FALSE,
              wtype = "multinom", simple = TRUE, drop0 = NA)
  })

  expect_equal(v0, vcov(fNA), tolerance = eps)

  set.seed(1234, "L")
  expect_no_condition({
    vNA <- vcovFWB(fit, R = 100, drop0 = TRUE,
                  wtype = "multinom")
  })

  expect_equal(vNA, vcov(fNA), tolerance = eps)

  # set.seed(1234, "L")
  # expect_no_condition({
  #   vs <- sandwich::vcovBS(fit, R = 100)
  # })
  #
  # expect_equal(v0, vs, tolerance = eps)

})
