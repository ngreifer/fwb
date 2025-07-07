test_that("p-values correctly invert CIs", {
  eps <- if (capabilities("long.double")) 1e-8 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  # Subsample to make BCa intervals faster
  set.seed(12345, "L")
  test_data <- test_data[sample.int(nrow(test_data), 500),]

  boot_fun <- function(data, w = NULL) {
    fit <- glm(Y_B ~ A + X1 + X2 + X3 + X4, data = data,
               family = quasibinomial("logit"), weights = w)
    coef(fit)
  }

  set.seed(1234, "L")
  expect_no_condition({
    boot <- fwb(test_data, boot_fun, R = 550, verbose = FALSE)
  })

  level <- .877
  index <- 4

  check_p_value_okay(boot, ci.type = "perc", level = level, index = index,
                     simultaneous = FALSE, eps = eps)

  check_p_value_okay(boot, ci.type = "bc", level = level, index = index,
                     simultaneous = FALSE, eps = eps)

  check_p_value_okay(boot, ci.type = "wald", level = level, index = index,
                     simultaneous = FALSE, eps = eps)

  check_p_value_okay(boot, ci.type = "norm", level = level, index = index,
                     simultaneous = FALSE, eps = eps)

  check_p_value_okay(boot, ci.type = "basic", level = level, index = index,
                     simultaneous = FALSE, eps = eps)

  check_p_value_okay(boot, ci.type = "bca", level = level, index = index,
                     simultaneous = FALSE, eps = eps)

  level <- .6
  index <- 2

  check_p_value_okay(boot, ci.type = "perc", level = level, index = index,
                     simultaneous = FALSE, eps = eps)

  check_p_value_okay(boot, ci.type = "bc", level = level, index = index,
                     simultaneous = FALSE, eps = eps)

  check_p_value_okay(boot, ci.type = "wald", level = level, index = index,
                     simultaneous = FALSE, eps = eps)

  check_p_value_okay(boot, ci.type = "norm", level = level, index = index,
                     simultaneous = FALSE, eps = eps)

  check_p_value_okay(boot, ci.type = "basic", level = level, index = index,
                     simultaneous = FALSE, eps = eps)

  check_p_value_okay(boot, ci.type = "bca", level = level, index = index,
                     simultaneous = FALSE, eps = eps)

})

test_that("p-values correctly invert simultaneous CIs", {
  skip_if_not_installed("mvtnorm")

  eps <- if (capabilities("long.double")) 1e-8 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  # Subsample to make BCa intervals faster
  set.seed(12345, "L")
  test_data <- test_data[sample.int(nrow(test_data), 500),]

  boot_fun <- function(data, w = NULL) {
    fit <- glm(Y_B ~ A + X1 + X2 + X3 + X4, data = data,
               family = quasibinomial("logit"), weights = w)
    coef(fit)
  }

  set.seed(1234, "L")
  expect_no_condition({
    boot <- fwb(test_data, boot_fun, R = 550, verbose = FALSE)
  })

  level <- .877
  index <- 4

  # check_p_value_okay(boot, ci.type = "perc", level = level, index = index,
  #                    simultaneous = TRUE, eps = eps)

  check_p_value_okay(boot, ci.type = "wald", level = level, index = index,
                     simultaneous = TRUE, eps = eps)

  level <- .6
  index <- 2

  # check_p_value_okay(boot, ci.type = "perc", level = level, index = index,
  #                    simultaneous = TRUE, eps = eps)

  check_p_value_okay(boot, ci.type = "wald", level = level, index = index,
                     simultaneous = TRUE, eps = eps)

})
