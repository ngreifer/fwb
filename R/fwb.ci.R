#' Fractional Weighted Bootstrap Confidence Intervals
#'
#' @description `fwb.ci()` generates several types of equi-tailed two-sided nonparametric confidence intervals. These include the normal approximation, the basic bootstrap interval, the percentile bootstrap interval, the bias-corrected percentile bootstrap interval, and the bias-correct and accelerated (BCa) bootstrap interval.
#'
#' @param fwb.out an `fwb` object; the output of a call to [fwb()].
#' @param conf the desired confidence level. Default is .95 for 95% confidence intervals.
#' @param type the type of confidence interval desired. Allowable options include `"norm"` (normal approximation), `"basic"` (basic interval), `"perc"` (percentile interval), `"bc"` (bias-correct percentile interval), and `"bca"` (BCa interval). More than one is allowed. Can also be `"all"` to request all of them. BCa intervals require that the number of bootstrap replications is larger than the sample size.
#' @param index the index of the position of the quantity of interest in `fwb.out$t0` if more than one was specified in `fwb()`. Only one value is allowed at a time. By default the first statistic is used.
#' @param h a function defining a transformation. The intervals are calculated on the scale of `h(t)` and the inverse function `hinv` applied to the resulting intervals. It must be a function of one variable only and for a vector argument, it must return a vector of the same length. Default is the identity function.
#' @param hinv a function, like `h`, which returns the inverse of `h`. It is used to transform the intervals calculated on the scale of `h(t)` back to the original scale. The default is the identity function. If `h` is supplied but `hinv` is not, then the intervals returned will be on the transformed scale.
#' @param ... ignored
#'
#' @returns
#' An `fwbci` object, which inherits from `bootci` and has the following components:
#' \item{R}{the number of bootstrap replications in the original call to `fwb()`.}
#' \item{t0}{the observed value of the statistic on the same scale as the intervals (i.e., after applying `h` and then `hinv`.}
#' \item{call}{the call to `fwb.ci()`.}
#'
#' There will be additional components named after each confidence interval type requested. For `"norm"`, this is a matrix with one row containing the confidence level and the two confidence interval limits. For the others, this is a matrix with one row containing the confidence level, the indices of the two order statistics used in the calculations, and the confidence interval limits.
#'
#' @details `fwb.ci()` functions similarly to \pkgfun{boot}{boot.ci} in that it takes in a bootstrapped object and computes confidence intervals. This interface is a bit old-fashioned, but was designed to mimic that of `boot.ci()`. For a more modern interface, see [summary.fwb()].
#'
#' The bootstrap intervals are defined as follows, with \eqn{\alpha =} 1 - `conf`, \eqn{t_0} the estimate in the original sample, \eqn{\hat{t}} the average of the bootstrap estimates, \eqn{s_t} the standard deviation of the bootstrap estimates, \eqn{t^{(i)}} the set of ordered estimates with \eqn{i} corresponding to their quantile, and \eqn{z_\frac{\alpha}{2}} and \eqn{z_{1-\frac{\alpha}{2}}} the upper and lower critical \eqn{z} scores.
#'
#' * `"norm"` (normal approximation): \eqn{[2t_0 - \hat{t} + s_t z_\frac{\alpha}{2}, 2t_0 - \hat{t} + s_t z_{1-\frac{\alpha}{2}}]}
#'
#' This involves subtracting the "bias" (\eqn{\hat{t} - t_0}) from the estimate \eqn{t_0} and using a standard Wald-type confidence interval. This method is valid when the statistic is normally distributed.
#'
#' * `"basic"`: \eqn{[2t_0 - t^{(1-\frac{\alpha}{2})}, 2t_0 - t^{(\frac{\alpha}{2})}]}
#'
#' * `"perc"` (percentile confidence interval): \eqn{[t^{(\frac{\alpha}{2})}, t^{(1-\frac{\alpha}{2})}]}
#'
#' * `"bc"` (bias-corrected percentile confidence interval): \eqn{[t^{(l)}, t^{(u)}]}
#'
#' \eqn{l = \Phi\left(2z_0 + z_\frac{\alpha}{2}\right)}, \eqn{u = \Phi\left(2z_0 + z_{1-\frac{\alpha}{2}}\right)}, where \eqn{\Phi(.)} is the normal cumulative density function (i.e., [pnorm()]) and \eqn{z_0 = \Phi^{-1}(q)} where \eqn{q} is the proportion of bootstrap estimates less than the original estimate \eqn{t_0}. This is similar to the percentile confidence interval but changes the specific quantiles of the bootstrap estimates to use, correcting for bias in the original estimate. It is described in Xu et al. (2020). When \eqn{t^0} is the median of the bootstrap distribution, the `"perc"` and `"bc"` intervals coincide.
#'
#' * `"bca"` (bias-corrected and accelerated confidence interval): \eqn{[t^{(l)}, t^{(u)}]}
#'
#' \eqn{l = \Phi\left(z_0 + \frac{z_0 + z_\frac{\alpha}{2}}{1-a(z_0+z_\frac{\alpha}{2})}\right)}, \eqn{u = \Phi\left(z_0 + \frac{z_0 + z_{1-\frac{\alpha}{2}}}{1-a(z_0+z_{1-\frac{\alpha}{2}})}\right)}, using the same definitions as above, but with the additional acceleration parameter \eqn{a}, where \eqn{a = \frac{1}{6}\frac{\sum{L^3}}{(\sum{L^2})^{3/2}}}. \eqn{L} is the empirical influence value of each unit, which is computed using the regression method described in \pkgfun{boot}{empinf}. When \eqn{a=0}, the `"bca"` and `"bc"` intervals coincide. The acceleration parameter corrects for bias and skewness in the statistic. It can only be used when clusters are absent and the number of bootstrap replications is larger than the sample size. Note that BCa intervals cannot be requested when `simple = TRUE` and there is randomness in thew `statistic` supplied to `fwb()`.
#'
#' Interpolation on the normal quantile scale is used when a non-integer order statistic is required, as in `boot::boot.ci()`. Note that unlike with `boot::boot.ci()`, studentized confidence intervals (`type = "stud"`) are not allowed.
#'
#' @seealso [fwb()] for performing the fractional weighted bootstrap; [get_ci()] for extracting confidence intervals from an `fwbci` object; [summary.fwb()] for producing clean output from `fwb()` that includes confidence intervals calculated by `fwb.ci()`; \pkgfun{boot}{boot.ci} for computing confidence intervals from the traditional bootstrap; [vcovFWB()] for computing parameter estimate covariance matrices using the fractional weighted bootstrap
#'
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
#' # Bias corrected percentile interval
#' bcci <- fwb.ci(fwb_out, index = "spontaneous", type = "bc")
#' bcci
#'
#' # Using `get_ci()` to extract confidence limits
#'
#' get_ci(bcci)
#'
#' # Interval calculated on original (log odds) scale,
#' # then transformed by exponentiation to be on OR
#' fwb.ci(fwb_out, index = "induced", type = "norm",
#'        hinv = exp)

#' @export
fwb.ci <- function(fwb.out, conf = .95, type = "bc", index = 1L,
                   h = base::identity, hinv = base::identity, ...) {

  call <- match.call()

  chk::chk_is(fwb.out, "boot")
  if (!chk::vld_number(conf) || conf <= .5 || conf >= 1) {
    .err("`conf` must be a single number between .5 and 1")
  }
  chk::chk_character(type)
  type <- match.arg(type, c("perc", "bc", "norm", "basic", "bca", "all"), several.ok = TRUE)

  if (any(type == "all")) {
    type <- c("perc", "bc", "norm", "basic", "bca")
  }

  if (any(type == "bca") &&
      fwb.out[["R"]] <= nrow(fwb.out[["data"]])) {
    msg <- "BCa confidence intervals cannot be computed when there are fewer bootstrap replications than units in the original dataset"

    if (all(type == "bca")) {
      .err(msg)
    }

    .wrn(msg)
    type <- type[type != "bca"]
  }

  if (any(type == "bca") &&
      isTRUE(attr(fwb.out, "simple", TRUE)) &&
      isTRUE(attr(fwb.out, "random_statistic", TRUE))) {
    msg <- "BCa confidence intervals cannot be computed when there is randomness in `statistic` and `simple = TRUE` in the call to `fbw()`. See `vignette(\"fwb-rep\")` for details"

    if (all(type == "bca")) {
      .err(msg)
    }

    .wrn(msg)
    type <- type[type != "bca"]
  }

  index <- check_index(index, fwb.out[["t"]])

  t0 <- fwb.out[["t0"]][index]
  t <- fwb.out[["t"]][, index]

  if (anyNA(t)) {
    .err("some bootstrap estimates are NA; cannot calculate confidence intervals")
  }

  if (!all(is.finite(t))) {
    .err("some bootstrap estimates are non-finite; cannot calculate confidence intervals")
  }

  if (all_the_same(t)) {
    .err(sprintf("all estimates are equal to %s\n Cannot calculate confidence intervals",
                 mean(t, na.rm = TRUE)))
  }

  if (length(t) != fwb.out[["R"]]) {
    .err(gettextf("'t' must be of length %d", fwb.out[["R"]]), domain = NA)
  }

  fins <- which(is.finite(t))
  t <- t[fins]
  R <- length(t)
  t0 <- h(t0)
  t <- h(t)

  output <- list(R = R, t0 = hinv(t0), call = call)

  if (any(type == "bc")) {
    output$bc <- bc.ci(t, t0, conf, hinv = hinv)
  }
  if (any(type == "perc")) {
    output$perc <- perc.ci(t, t0, conf, hinv = hinv)
  }
  if (any(type == "norm")) {
    output$norm <- norm.ci(t, t0, conf, hinv = hinv)
  }
  if (any(type == "basic")) {
    output$basic <- basic.ci(t, t0, conf, hinv = hinv)
  }
  if (any(type == "bca")) {
    output$bca <- bca.ci(t, t0, fwb.out, index, conf, hinv = hinv, h = h)
  }

  attr(output, "conf") <- conf

  class(output) <- c("fwbci", "bootci")

  output
}

#' @describeIn fwb.ci Print a bootstrap confidence interval
#'
#' @param x an `fwbci` object; the output of a call to `fwb.ci()`.
#'
#' @exportS3Method print fwbci
print.fwbci <- function (x, hinv = NULL, ...) {
  ci.out <- x
  cl <- ci.out[["call"]]
  ntypes <- length(ci.out) - 3L
  nints <- nrow(ci.out[[4L]])
  t0 <- ci.out[["t0"]]
  if (!is.null(hinv))
    t0 <- hinv(t0)
  digs <- ceiling(log10(abs(t0)))
  if (digs <= 0)
    digs <- 4
  else if (digs >= 4)
    digs <- 0
  else digs <- 4 - digs
  intlabs <- NULL
  basrg <- bcperg <- perg <- bcarg <- NULL
  if (!is.null(ci.out[["norm"]])) {
    intlabs <- c(intlabs, "      Normal   ")
  }
  if (!is.null(ci.out[["basic"]])) {
    intlabs <- c(intlabs, "       Basic   ")
    basrg <- range(ci.out[["basic"]][, 2:3])
  }
  if (!is.null(ci.out[["bc"]])) {
    intlabs <- c(intlabs, "BC Percentile  ")
    bcperg <- range(ci.out[["bc"]][, 2:3])
  }
  if (!is.null(ci.out[["perc"]])) {
    intlabs <- c(intlabs, "  Percentile   ")
    perg <- range(ci.out[["perc"]][, 2:3])
  }
  if (!is.null(ci.out[["bca"]])) {
    intlabs <- c(intlabs, "        BCa    ")
    bcarg <- range(ci.out[["bca"]][, 2:3])
  }

  level <- 100 * attr(ci.out, "conf", TRUE)
  if (ntypes == 4L)
    n1 <- n2 <- 2L
  else if (ntypes == 5L) {
    n1 <- 3L
    n2 <- 2L
  }
  else {
    n1 <- ntypes
    n2 <- 0L
  }
  ints1 <- matrix(NA, nints, 2L * n1 + 1L)
  ints1[, 1L] <- level
  n0 <- 4L
  for (i in n0:(n0 + n1 - 1)) {
    j <- c(2L * i - 6L, 2L * i - 5L)
    nc <- ncol(ci.out[[i]])
    nc <- c(nc - 1L, nc)
    if (is.null(hinv)) ints1[, j] <- ci.out[[i]][, nc]
    else ints1[, j] <- hinv(ci.out[[i]][, nc])
  }
  n0 <- 4L + n1
  ints1 <- format(round(ints1, digs))
  ints1[, 1L] <- paste("\n", level, "%  ", sep = "")
  ints1[, 2 * (1L:n1)] <- paste("(", ints1[, 2 * (1L:n1)], ",", sep = "")
  ints1[, 2 * (1L:n1) + 1L] <- paste(ints1[, 2 * (1L:n1) + 1L], ")  ")
  if (n2 > 0) {
    ints2 <- matrix(NA, nints, 2L * n2 + 1L)
    ints2[, 1L] <- level
    j <- c(2L, 3L)
    for (i in n0:(n0 + n2 - 1L)) {
      nc <- ncol(ci.out[[i]])
      nc <- c(nc - 1L, nc)
      if (is.null(hinv))
        ints2[, j] <- ci.out[[i]][, nc]
      else ints2[, j] <- hinv(ci.out[[i]][, nc])
      j <- j + 2L
    }
    ints2 <- format(round(ints2, digs))
    ints2[, 1L] <- paste("\n", level, "%  ", sep = "")
    ints2[, 2 * (1L:n2)] <- paste("(", ints2[, 2 * (1L:n2)], ",", sep = "")
    ints2[, 2 * (1L:n2) + 1L] <- paste(ints2[, 2 * (1L:n2) + 1L], ")  ")
  }
  R <- ci.out[["R"]]
  cat("BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS\n")
  cat(paste("Based on", R, "bootstrap replicates\n\n"))
  cat("CALL : \n")
  dput(cl, control = NULL)
  cat("\nIntervals : ")
  cat("\nLevel", intlabs[1L:n1])
  cat(t(ints1))
  if (n2 > 0) {
    cat("\n\nLevel", intlabs[(n1 + 1):(n1 + n2)])
    cat(t(ints2))
  }
  if (!is.null(cl[["h"]])) {
    if (is.null(cl[["hinv"]]) && is.null(hinv))
      cat("\nCalculations and Intervals on ", "Transformed Scale\n")
    else cat("\nCalculations on Transformed Scale;", " Intervals on Original Scale\n")
  }
  else if (is.null(cl[["hinv"]]) && is.null(hinv))
    cat("\nCalculations and Intervals on Original Scale\n")
  else cat("\nCalculations on Original Scale", " but Intervals Transformed\n")

  if (!is.null(basrg)) {
    if ((basrg[1L] <= 1) || (basrg[2L] >= R))
      cat("Warning : Basic Intervals used Extreme Quantiles\n")
    if ((basrg[1L] <= 10) || (basrg[2L] >= R - 9))
      cat("Some basic intervals may be unstable\n")
  }
  if (!is.null(bcperg)) {
    if ((bcperg[1L] <= 1) || (bcperg[2L] >= R))
      cat("Warning : Bias-Corrected Percentile Intervals used Extreme Quantiles\n")
    if ((bcperg[1L] <= 10) || (bcperg[2L] >= R - 9))
      cat("Some bias-corrected percentile intervals may be unstable\n")
  }
  if (!is.null(perg)) {
    if ((perg[1L] <= 1) || (perg[2L] >= R))
      cat("Warning : Percentile Intervals used Extreme Quantiles\n")
    if ((perg[1L] <= 10) || (perg[2L] >= R - 9))
      cat("Some percentile intervals may be unstable\n")
  }
  if (!is.null(bcarg)) {
    if ((bcarg[1L] <= 1) || (bcarg[2L] >= R))
      cat("Warning : BCa Intervals used Extreme Quantiles\n")
    if ((bcarg[1L] <= 10) || (bcarg[2L] >= R - 9))
      cat("Some BCa intervals may be unstable\n")
  }

  invisible(ci.out)
}

#' @title Extract Confidence Intervals from a `bootci` Object
#'
#' @description `get_ci()` extracts the confidence intervals from the output of a call to \pkgfun{boot}{boot.ci} or [fwb.ci()] in a clean way. Normally the confidence intervals can be a bit challenging to extract because of the unusual structure of the object.
#'
#' @param x a `bootci` object; the output of a call to `boot::boot.ci()` or `fwb.ci()`.
#' @param type the type of confidence intervals to extract. Only those available in `x` are allowed. Should be a given as a subset of the types passed to `type` in `boot.ci()` or `fwb.ci()`. The default, `"all"`, extracts all confidence intervals in `x`.
#'
#' @returns
#' A list with an entry for each confidence interval type; each entry is a numeric vector of length 2 with names `"L"` and `"U"` for the lower and upper interval bounds, respectively. The `"conf"` attribute contains the confidence level.
#'
#' @examples
#' #See example at help("fwb.ci")
#'
#' @seealso [fwb.ci()], [confint.fwb()], \pkgfun{boot}{boot.ci}

#' @export
get_ci <- function(x, type = "all") {
  chk::chk_is(x, "bootci")
  chk::chk_character(type)

  name_trans <- function(x, old, new) {
    if (old %in% names(x)) names(x)[names(x) == old] <- new
    x
  }

  x <- name_trans(x, "normal", "norm")
  x <- name_trans(x, "student", "stud")
  x <- name_trans(x, "percent", "perc")

  allowed_cis <- intersect(c("perc", "bc", "stud", "norm", "basic", "bca"),
                           names(x))

  if ("all" %in% type) {
    type <- allowed_cis
  }
  else {
    type <- match_arg(type, allowed_cis)
  }

  out <- setNames(lapply(type, function(t) {
    nc <- ncol(x[[t]])
    setNames(x[[t]][, c(nc - 1, nc)], c("L", "U"))
  }), type)

  attr(out, "conf") <- {
    if (is_not_null(attr(x, "conf", TRUE))) attr(x, "conf", TRUE)
    else x[[4L]][, 1L]
  }

  out
}

perc.ci <- function(t, t0, conf = 0.95, hinv = identity) {
  alpha <- (1 + c(-conf, conf))/2
  qq <- norm.inter(t, alpha)
  cbind(conf, matrix(qq[, 1L], ncol = 2L), matrix(hinv(qq[, 2]), ncol = 2L))
}

norm.inter <- function(t, alpha) {
  t <- t[is.finite(t)]
  R <- length(t)
  rk <- (R + 1) * alpha
  if (!all(rk > 1 & rk < R)) {
    .wrn("extreme order statistics used as endpoints")
  }
  k <- trunc(rk)
  inds <- seq_along(k)
  out <- inds
  kvs <- k[k > 0 & k < R]
  tstar <- sort(t, partial = sort(union(c(1, R), c(kvs, kvs + 1))))
  ints <- (k == rk)
  if (any(ints)) {
    out[inds[ints]] <- tstar[k[inds[ints]]]
  }
  out[k == 0] <- tstar[1L]
  out[k == R] <- tstar[R]
  not <- function(v) xor(rep(TRUE, length(v)), v)
  temp <- inds[not(ints) & k != 0 & k != R]
  temp1 <- qnorm(alpha[temp])
  temp2 <- qnorm(k[temp]/(R + 1))
  temp3 <- qnorm((k[temp] + 1)/(R + 1))
  tk <- tstar[k[temp]]
  tk1 <- tstar[k[temp] + 1L]
  out[temp] <- tk + (temp1 - temp2)/(temp3 - temp2) * (tk1 - tk)
  cbind(round(rk, 2), out)
}

bc.ci <- function(t, t0, conf = .95, hinv = identity) {
  alpha <- (1 + c(-conf, conf))/2
  s <- min(max(1, sum(t < t0)), length(t) - 1)
  w <- qnorm(s/length(t))

  zalpha <- qnorm(alpha)

  adj.alpha <- pnorm(w + (w + zalpha))

  qq <- norm.inter(t, adj.alpha)
  cbind(conf, matrix(qq[, 1L], ncol = 2L), matrix(hinv(qq[, 2]), ncol = 2L))
}

norm.ci <- function(t, t0, conf = 0.95, hinv = identity) {
  bias <- mean(t) - t0

  merr <- sd(t) * qnorm((1 + conf)/2)
  cbind(conf, hinv(t0 - bias - merr), hinv(t0 - bias + merr))
}

basic.ci <- function (t, t0, conf = 0.95, hinv = identity) {
  alpha <- (1 + c(-conf, conf))/2
  qq <- norm.inter(t, rev(alpha))
  cbind(conf,
        matrix(qq[, 1L], ncol = 2L),
        matrix(hinv(2 * t0 - qq[, 2L]), ncol = 2L))
}

bca.ci <- function(t, t0, boot.out, index, conf = .95, hinv = identity, h = identity) {
  if (is_not_null(boot.out[["cluster"]])) {
    .err("the BCa confidence interval cannot be used with clusters")
  }

  alpha <- (1 + c(-conf, conf))/2
  s <- min(max(1, sum(t < t0)), length(t) - 1)
  w <- qnorm(s/length(t))

  zalpha <- qnorm(alpha)

  L <- empinf.reg(boot.out, t = t)

  a <- sum(L^3)/(6 * sum(L^2)^1.5)

  if (!is.finite(a)) {
    .err("estimated adjustment 'a' is NA")
  }

  adj.alpha <- pnorm(w + (w + zalpha)/(1 - a * (w + zalpha)))
  qq <- norm.inter(t, adj.alpha)

  cbind(conf, matrix(qq[, 1L], ncol = 2L), matrix(hinv(h(qq[, 2L])), ncol = 2L))
}

empinf.reg <- function(boot.out, t) {
  fins <- which(is.finite(t))
  t <- t[fins]
  n <- NROW(boot.out[["data"]])
  f <- boot.array(boot.out)[fins, ]
  X <- f/n
  X[, 1L] <- 1
  beta <- .lm.fit(x = X, y = t)$coefficients[-1L]
  l <- rep.int(0, n)
  l[-1L] <- beta

  l - mean(l)
}

boot.array <- function(boot.out) {
  if (identical(attr(boot.out, "boot_type", TRUE), "boot")) {
    rlang::check_installed("boot")
    return(boot::boot.array(boot.out))
  }

  genv <- globalenv()

  #Return seed to its prior state after generating weights using seed from boot.out
  old_seed <- get(".Random.seed", envir = genv, inherits = FALSE)
  on.exit(suspendInterrupts({
    if (is_null(old_seed)) {
      rm(".Random.seed", envir = genv, inherits = FALSE)
    }
    else {
      assign(".Random.seed", value = old_seed, envir = genv, inherits = FALSE)
    }
  }))

  assign(".Random.seed", value = boot.out[["seed"]], envir = genv)

  gen_weights <- make_gen_weights(boot.out[["wtype"]])

  n <- nrow(boot.out[["data"]])
  R <- boot.out[["R"]]

  if (!isTRUE(attr(boot.out, "simple", TRUE)) || is_null(attr(boot.out, "cl", TRUE))) {
    return(gen_weights(n, R, boot.out[["strata"]]))
  }

  if (isTRUE(attr(boot.out, "simple", TRUE)) &&
      isTRUE(attr(boot.out, "random_statistic", TRUE))) {
    .wrn("bootstrap weights cannot be repliably re-generated when there is randomness in `statistic` and `simple = TRUE` in the call to `fbw()`. See `vignette(\"fwb-rep\")` for details")
  }

  FUN <- function(i) {
    drop(gen_weights(n, 1L, boot.out[["strata"]]))
  }

  opb <- pbapply::pboptions(type = "none")
  on.exit(pbapply::pboptions(opb))

  #Run bootstrap
  if (identical(attr(boot.out, "cl", TRUE), "future"))
    do.call("rbind", pbapply::pblapply(seq_len(R), FUN, cl = "future", future.seed = TRUE))
  else
    do.call("rbind", pbapply::pblapply(seq_len(R), FUN, cl = attr(boot.out, "cl", TRUE)))
}


