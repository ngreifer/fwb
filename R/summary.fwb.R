#' Summarize `fwb` Output
#'
#' `summary()` creates a regression summary-like table that displays the bootstrap estimates, their empirical standard errors, and their confidence intervals. `confint()` produces just the confidence intervals, which are computed using [fwb.ci()], and is called internally by `summary()`.
#'
#' @param object an `fwb` object; the output of a call to [fwb()].
#' @param conf,level the desired confidence level. Default is .95 for 95% confidence intervals.
#' @param ci.type the type of confidence interval desired. Allowable options include `"norm"` (normal approximation), `"basic"` (basic interval), `"perc"` (percentile interval), `"bc"` (bias-correct percentile interval), and `"bca"` (bias-corrected and accelerated \[BCa\] interval). Only one is allowed. BCa intervals require that the number of bootstrap replications is larger than the sample size. See [fwb.ci()] for details. The default is `"bc"`.
#' @param index,parm the index or indices of the position of the quantity of interest in `x$t0` if more than one was specified in `fwb()`. Default is to display all quantities.
#' @param p.value `logical`; whether to display p-values for the test that each parameter is equal to 0. The p-value is computed using a Z-test with the test statistic computed as the ratio of the estimate to its bootstrap standard error. This test is only valid when the bootstrap distribution is normally distributed around 0 and is not guaranteed to agree with any of the confidence intervals. Default is `FALSE`.
#' @param ... ignored.
#'
#' @return Fpr `summary()`, a `summary.fwb` object, which is a matrix with the following columns:
#' * `Estimate`: the statistic estimated in the original sample
#' * `Std. Error`: the standard deviation of the bootstrap estimates
#' * `CI {L}%` and `CI {U}%`, the upper and lower confidence interval bounds computed using the argument to `ci.type`.
#'
#' When `p.value = TRUE`, two additional columns, `z value` and `Pr(>|z|)`, are included containing the z-statistic and p-value for each computed statistic, respectively.
#'
#' For `confint()`, a matrix with a row for each statistic and a column for the upper and lower confidence interval limits.
#'
#' @seealso [fwb()] for performing the fractional weighted bootstrap; [fwb.ci()] for computing multiple confidence intervals for a single bootstrapped quantity
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
#' fwb_out <- fwb(infert, fit_fun, R = 199, verbose = FALSE)
#'
#' # Basic confidence interval for both estimates
#' summary(fwb_out, ci.type = "basic")
#'
#' # Just for "induced" coefficient; p-values requested
#' summary(fwb_out, index = "induced", p.value = TRUE)
#'

#' @exportS3Method summary fwb
summary.fwb <- function(object, conf = .95, ci.type = "bc", p.value = FALSE, index = 1L:ncol(object$t), ...) {

  chk::chk_number(conf)
  chk::chk_lt(conf, 1)
  chk::chk_gt(conf, .5)
  chk::chk_string(ci.type)
  chk::chk_flag(p.value)

  index <- check_index(index, object[["t"]], several.ok = TRUE)

  ci.type <- match_arg(ci.type, c("perc", "bc", "norm", "basic", "bca"))

  pct <- fmt.prc(c((1 - conf) / 2, 1 - (1 - conf) / 2))

  out <- matrix(nrow = length(index), ncol = 4L,
                dimnames = list(names(object[["t0"]])[index],
                                c("Estimate", "Std. Error", paste("CI", pct))))

  out[, "Estimate"] <- object[["t0"]][index]
  out[, "Std. Error"] <- apply(object[["t"]][, index, drop = FALSE], 2L, sd)

  out[, 3:4] <- confint.fwb(object, parm = index, level = conf, ci.type = ci.type)

  if (p.value) {
    z <- out[, "Estimate"] / out[, "Std. Error"]
    out <- cbind(out, `z value` = z, `Pr(>|z|)` = 2 * pnorm(-abs(z)))
  }

  class(out) <- c("summary.fwb", class(out))
  out
}

#' @exportS3Method confint fwb
#' @rdname summary.fwb
confint.fwb <- function(object, parm, level = .95, ci.type = "bc", ...) {

  chk::chk_number(level)
  chk::chk_lt(level, 1)
  chk::chk_gt(level, .5)
  chk::chk_string(ci.type)

  if (missing(parm)) {
    parm <- seq_len(ncol(object$t))
  }

  index <- check_index(parm, object[["t"]], several.ok = TRUE)

  ci.type <- match_arg(ci.type, c("perc", "bc", "norm", "basic", "bca"))

  ci.list <- lapply(index, function(i) {
    fwb.ci(object, conf = level, type = ci.type, index = i)
  })

  pct <- fmt.prc(c((1 - level) / 2, 1 - (1 - level) / 2))

  out <- matrix(NA_real_, nrow = length(index), ncol = 2L,
                dimnames = list(names(object[["t0"]])[index],
                                pct))

  out[, 1L] <- vapply(ci.list, function(x) {
    nc <- ncol(x[[ci.type]])
    x[[ci.type]][nc - 1L]
  }, numeric(1L))

  out[, 2L] <- vapply(ci.list, function(x) {
    nc <- ncol(x[[ci.type]])
    x[[ci.type]][nc]
  }, numeric(1L))

  out
}

#' @exportS3Method print summary.fwb
print.summary.fwb <- function(x, digits = 3L, ...) {
  has.p <- ncol(x) > 4L
  stats::printCoefmat(x, digits = digits, cs.ind = 1:4,
                      tst.ind = if (has.p) 5L else integer(),
                      has.Pvalue = has.p)
  invisible(x)
}
