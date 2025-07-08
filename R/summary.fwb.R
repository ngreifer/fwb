#' Summarize `fwb` Output
#'
#' `summary()` creates a regression summary-like table that displays the bootstrap estimates, their empirical standard errors, their confidence intervals, and, optionally, p-values for tests against a null value. `confint()` produces just the confidence intervals, and is called internally by `summary()`.
#'
#' @param object an `fwb` object; the output of a call to [fwb()].
#' @param conf,level the desired confidence level. Default is .95 for 95% confidence intervals. Set to 0 to prevent calculation of confidence intervals.
#' @param ci.type the type of confidence interval desired. Allowable options include `"wald"` (Wald interval), `"norm"` (normal approximation), `"basic"` (basic interval), `"perc"` (percentile interval), `"bc"` (bias-corrected percentile interval), and `"bca"` (bias-corrected and accelerated \[BCa\] interval). Only one is allowed. BCa intervals require the number of bootstrap replications to be larger than the sample size. See [fwb.ci()] for details. The default is `"bc"`. Ignored if both `conf = 0` and `p.values = FALSE`.
#' @param index,parm the index or indices of the position of the quantity of interest if more than one was specified in `fwb()`. Default is to display all quantities.
#' @param p.value `logical`; whether to display p-values for the test that each parameter is equal to `null`. Default is `FALSE`. See Details.
#' @param null `numeric`; when `p.value = TRUE`, the value of the estimate under the null hypothesis. Default is 0. Only one value can be supplied and it is applied to all tests.
#' @param simultaneous `logical`; whether to adjust confidence intervals and p-values to ensure correct familywise coverage/size across all specified estimates. See Details. Default is `FALSE` for standard pointwise intervals. `TRUE` is only allowed when `ci.type` is `"wald"` or `"perc"`.
#' @param ... ignored.
#'
#' @return
#' For `summary()`, a `summary.fwb` object, which is a matrix with the following columns:
#' * `Estimate`: the statistic estimated in the original sample
#' * `Std. Error`: the standard deviation of the bootstrap estimates
#' * `CI {L}%` and `CI {U}%`: the upper and lower confidence interval bounds computed using the argument to `ci.type` (only when `conf` is not 0).
#' * `z value`: when `p.value = TRUE` and `ci.type = "wald"`, the z-statistic for the test of the estimate against against `null`.
#' * `Pr(>|z|)`: when `p.value = TRUE`, the p-value for the test of the estimate against against `null`.
#'
#' For `confint()`, a matrix with a row for each statistic and a column for the upper and lower confidence interval limits.
#'
#' @details
#' P-values are computed by inverting the confidence interval for each parameter, i.e., finding the largest confidence level yielding a confidence interval that excludes `null`, and taking the p-value to be one minus that level. This ensures conclusions from tests based on the p-value and whether the confidence interval contains the null value always yield the same conclusion. Prior to version 0.5.0, all p-values were based on inverting Wald confidence intervals, regardless of `ci.type`.
#'
#' Simultaneous confidence intervals are computed using the "sup-t" confidence band, which involves modifying the confidence level so that the intersection of all the adjusted confidence intervals contain the whole parameter vector with the specified coverage. This will always be less conservative than Bonferroni or Holm adjustment. See Olea and Plagborg-Møller (2019) for details on implementation for Wald and percentile intervals. Simultaneous p-values are computed by inverting the simultaneous bands. Simultaneous inference is only allowed when `ci.type` is `"wald"` or `"perc"` and `index` has length greater than 1. When `ci.type = "wald"`, the \pkg{mvtnorm} package must be installed.
#'
#' `tidy()` and `print()` methods are available for `summary.fwb` objects.
#'
#' @seealso [fwb()] for performing the fractional weighted bootstrap; [fwb.ci()] for computing multiple confidence intervals for a single bootstrapped quantity
#'
#' @references
#' Montiel Olea, J. L., & Plagborg-Møller, M. (2019). Simultaneous confidence bands: Theory, implementation, and an application to SVARs. *Journal of Applied Econometrics*, 34(1), 1–17. \doi{10.1002/jae.2656}
#'
#' @examples
#' set.seed(123, "L'Ecuyer-CMRG")
#' data("infert")
#'
#' fit_fun <- function(data, w) {
#'   fit <- glm(case ~ spontaneous + induced, data = data,
#'              family = "quasibinomial", weights = w)
#'   coef(fit)
#' }
#'
#' fwb_out <- fwb(infert, fit_fun, R = 199,
#'                verbose = FALSE)
#'
#' # Basic confidence interval for both estimates
#' summary(fwb_out, ci.type = "basic")
#'
#' # Just for "induced" coefficient; p-values requested,
#' # no confidence intervals
#' summary(fwb_out, ci.type = "norm", conf = 0,
#'         index = "induced", p.value = TRUE)
#'

#' @exportS3Method summary fwb
summary.fwb <- function(object, conf = .95, ci.type = "bc", p.value = FALSE,
                        index = seq_len(ncol(object$t)), null = 0, simultaneous = FALSE, ...) {

  chk::chk_number(conf)
  chk::chk_range(conf, c(0, 1), inclusive = TRUE)

  chk::chk_flag(p.value)

  index <- check_index(index, object[["t"]], several.ok = TRUE)

  if (p.value || conf > 0) {
    chk::chk_string(ci.type)
    ci.type <- match_arg(ci.type, .allowed_ci.types())

    if (length(index) <= 1L) {
      simultaneous <- FALSE
    }

    chk::chk_flag(simultaneous)

    if (simultaneous) {
      if (!ci.type %in% c("perc", "wald")) {
      .err('simultaneous inference can only be used when `ci.type` is `"wald"` or `"perc"`')
      }

      if (ci.type == "wald" && conf > 0 && conf <= .5) {
        .err('`conf` must be greater than .5 to use with `simultaneous = TRUE` and `ci.type = "wald"`')
      }
    }
  }

  out <- matrix(nrow = length(index), ncol = 2L,
                dimnames = list(names(object[["t0"]])[index],
                                c("Estimate", "Std. Error")))

  out[, "Estimate"] <- object[["t0"]][index]
  out[, "Std. Error"] <- apply(object[["t"]][, index, drop = FALSE], 2L, sd)

  if (conf > 0) {
    chk::chk_lt(conf, 1)

    pct <- fmt.prc(c((1 - conf) / 2, 1 - (1 - conf) / 2))

    ci <- confint.fwb(object, parm = index, level = conf, ci.type = ci.type,
                      simultaneous = simultaneous)

    colnames(ci) <- paste("CI", pct)

    out <- cbind(out, ci)
  }

  if (p.value) {
    chk::chk_number(null)

    if (ci.type == "wald") {
      z <- (out[, "Estimate"] - null) / out[, "Std. Error"]
      out <- cbind(out, `z value` = z)
    }

    p.values <- invert_ci(ci.type, t = object[["t"]], t0 = object[["t0"]],
                          null = null, index = index, boot.out = object)

    if (simultaneous) {
      p.values <- simultaneous_p_value(object, p.values, index, ci.type)
    }

    out <- cbind(out, `Pr(>|z|)` = p.values)
  }

  attr(out, "conf") <- conf

  if (p.value || conf > 0) {
    attr(out, "ci.type") <- ci.type
    attr(out, "simultaneous") <- simultaneous
  }

  if (p.value) {
    attr(out, "null") <- null
  }

  class(out) <- c("summary.fwb", class(out))

  out
}

#' @exportS3Method confint fwb
#' @rdname summary.fwb
confint.fwb <- function(object, parm, level = .95, ci.type = "bc", simultaneous = FALSE, ...) {

  chk::chk_number(level)
  chk::chk_range(level, c(0, 1), inclusive = FALSE)
  chk::chk_string(ci.type)

  if (missing(parm)) {
    parm <- seq_len(ncol(object$t))
  }

  index <- check_index(parm, object[["t"]], several.ok = TRUE)

  if (length(index) <= 1L) {
    simultaneous <- FALSE
  }

  chk::chk_flag(simultaneous)

  ci.type <- match_arg(ci.type, .allowed_ci.types())

  if (simultaneous) {
    if (!ci.type %in% c("perc", "wald")) {
      .err('simultaneous inference can only be used when `ci.type` is `"wald"` or `"perc"`')
    }

    if (ci.type == "wald" && level <= .5) {
      .err('`level` must be greater than .5 to use with `simultaneous = TRUE` and `ci.type = "wald"`')
    }

    new_level <- simultaneous_ci_level(object, level, index, ci.type)
  }
  else {
    new_level <- level
  }

  ci.out <- compute_ci(ci.type, t = object[["t"]], t0 = object[["t0"]],
                       conf = new_level, index = index, boot.out = object)

  pct <- fmt.prc(c((1 - level) / 2, 1 - (1 - level) / 2))

  out <- matrix(NA_real_, nrow = length(index), ncol = 2L,
                dimnames = list(names(object[["t0"]])[index],
                                pct))

  nc <- ncol(ci.out)

  out[, 1:2] <- ci.out[, c(nc - 1L, nc)]

  attr(out, "simultaneous") <- simultaneous

  out
}

#' @exportS3Method print summary.fwb
print.summary.fwb <- function(x, digits = 3L, ...) {
  has.p <- is_not_null(attr(x, "null", TRUE))
  has.ci <- attr(x, "conf", TRUE) > 0
  has.z <- has.p && identical(attr(x, "ci.type", TRUE), "wald")

  stats::printCoefmat(x, digits = digits,
                      cs.ind = if (has.ci) 1:4 else 1:2,
                      tst.ind = if (has.z && has.ci) 5L else if (has.z) 3L else integer(),
                      has.Pvalue = has.p,
                      na.print = ".")
  invisible(x)
}

#' @importFrom generics tidy
#' @exportS3Method tidy summary.fwb
tidy.summary.fwb <- function(x, ...) {
  has.p <- is_not_null(attr(x, "null", TRUE))
  has.ci <- attr(x, "conf", TRUE) > 0
  has.z <- has.p && identical(attr(x, "ci.type", TRUE), "wald")

  x <- as.data.frame(x)

  names(x)[1:2] <- c("estimate", "std.error")

  if (has.ci) {
    names(x)[3:4] <- c("conf.low", "conf.high")
  }

  if (has.z) {
    if (has.ci) {
      names(x)[5L] <- "statistic"
    }
    else {
      names(x)[3L] <- "statistic"
    }
  }

  if (has.p) {
    names(x)[ncol(x)] <- "p.value"
  }

  x <- cbind(rownames(x), x)
  names(x)[1L] <- "term"

  rownames(x) <- NULL

  class(x) <- c("tbl_df", "tbl", "data.frame")

  x
}
