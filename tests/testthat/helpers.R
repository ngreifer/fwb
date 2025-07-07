expect_not_equal <- function(object, expected, ...,
                             tolerance = if (edition_get() >= 3) testthat_tolerance(),
                             info = NULL, label = NULL, expected.label = NULL) {

  if (!capabilities("long.double")) {
    return(NULL)
  }

  act <- quasi_label(rlang::enquo(object), label, arg = "object")
  exp <- quasi_label(rlang::enquo(expected), expected.label, arg = "expected")

  if (edition_get() >= 3) {
    expect_waldo_not_equal("equal", act, exp, info, ..., tolerance = tolerance)
  }
  else {
    if (!is.null(tolerance)) {
      comp <- compare(act$val, exp$val, ..., tolerance = tolerance)
    }
    else {
      comp <- compare(act$val, exp$val, ...)
    }
    expect(!comp$equal, sprintf("%s equal to %s.\n%s", act$lab, exp$lab, comp$message), info = info)
    invisible(act$val)
  }
}

expect_waldo_not_equal <- function(type, act, exp, info, ...) {
  comp <- waldo::compare(act$val, exp$val, ..., x_arg = "actual",
                         y_arg = "expected")

  expect(length(comp) > 0, sprintf("%s (%s) is %s to %s (%s).\n\n%s",
                                   act$lab, "`actual`", type, exp$lab,
                                   "`expected`",
                                   paste0(comp, collapse = "\n\n")),
         info = info, trace_env = rlang::caller_env())
  invisible(act$val)
}

#Use regex to make strings invariant to white spaces
.w <- function(x) {
  gsub(" ", "(\\s+)", x, fixed = TRUE)
}

check_p_value_okay <- function(boot, ci.type, level, index, simultaneous = FALSE, eps) {
  R <- boot[["R"]]
  if (simultaneous) {

    eps_ <- max(eps, 1e-3)

    set.seed(7)
    suppressWarnings({
      s0 <- summary(boot, conf = level, ci.type = ci.type,
                   simultaneous = TRUE)
    })

    expect_s3_class(s0, "summary.fwb")
    expect_true(attr(s0, "simultaneous", TRUE))

    ## lower bound
    suppressWarnings({
      p <- summary(boot, conf = 0, p.value = TRUE, ci.type = ci.type,
                   simultaneous = TRUE,
                   null = s0[index, 3])[index, "Pr(>|z|)"]
    })

    suppressWarnings({
      s <- summary(boot, conf = 1 - p, ci.type = ci.type,
                   simultaneous = TRUE)
    })

    expect_equal(s[index, 3:4], s0[index, 3:4],
                 tolerance = eps_,
                 ignore_attr = TRUE)

    ## upper bound
    suppressWarnings({
      p <- summary(boot, conf = 0, p.value = TRUE, ci.type = ci.type,
                   simultaneous = TRUE,
                   null = s0[index, 4])[index, "Pr(>|z|)"]
    })

    suppressWarnings({
      s <- summary(boot, conf = 1 - p, ci.type = ci.type,
                   simultaneous = TRUE)
    })

    expect_equal(s[index, 3:4], s0[index, 3:4],
                 tolerance = eps_,
                 ignore_attr = TRUE)
  }
  else {
    suppressWarnings({
      s <- summary(boot, conf = level, ci.type = ci.type,
                   index = index)
    })

    expect_s3_class(s, "summary.fwb")
    expect_false(attr(s, "simultaneous", TRUE))

    ## lower bound
    suppressWarnings({
      p <- summary(boot, conf = 0, p.value = TRUE, ci.type = ci.type,
                   index = index, null = s[1, 3])[1, "Pr(>|z|)"]
    })

    expect_equal(p, 1 - level, tolerance = eps)

    ## upper bound
    suppressWarnings({
      p <- summary(boot, conf = 0, p.value = TRUE, ci.type = ci.type,
                   index = index, null = s[1, 4])[1, "Pr(>|z|)"]
    })

    expect_equal(p, 1 - level, tolerance = eps)
  }
}
