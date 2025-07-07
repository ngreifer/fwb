#' Fractional Weighted Bootstrap Confidence Intervals
#'
#' @description `fwb.ci()` generates several types of equi-tailed two-sided nonparametric confidence intervals. These include the normal approximation, the basic bootstrap interval, the percentile bootstrap interval, the bias-corrected percentile bootstrap interval, and the bias-correct and accelerated (BCa) bootstrap interval.
#'
#' @param fwb.out an `fwb` object; the output of a call to [fwb()].
#' @param conf the desired confidence level. Default is .95 for 95% confidence intervals.
#' @param type the type of confidence interval desired. Allowable options include `"wald"` (Wald interval), `"norm"` (normal approximation), `"basic"` (basic interval), `"perc"` (percentile interval), `"bc"` (bias-correct percentile interval), and `"bca"` (BCa interval). More than one is allowed. Can also be `"all"` to request all of them. BCa intervals require that the number of bootstrap replications is larger than the sample size.
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
#' There will be additional components named after each confidence interval type requested. For `"wald"` and `"norm"`, this is a matrix with one row containing the confidence level and the two confidence interval limits. For the others, this is a matrix with one row containing the confidence level, the indices of the two order statistics used in the calculations, and the confidence interval limits.
#'
#' @details `fwb.ci()` functions similarly to \pkgfun{boot}{boot.ci} in that it takes in a bootstrapped object and computes confidence intervals. This interface is a bit old-fashioned, but was designed to mimic that of `boot.ci()`. For a more modern interface, see [summary.fwb()].
#'
#' The bootstrap intervals are defined as follows, with \eqn{\alpha =} 1 - `conf`, \eqn{t_0} the estimate in the original sample, \eqn{\hat{t}} the average of the bootstrap estimates, \eqn{s_t} the standard deviation of the bootstrap estimates, \eqn{t^{(i)}} the set of ordered estimates with \eqn{i} corresponding to their quantile, and \eqn{z_\frac{\alpha}{2}} and \eqn{z_{1-\frac{\alpha}{2}}} the upper and lower critical \eqn{z} scores.
#'
#' \describe{
#' \item{`"wald"`}{
#'   \deqn{\left[t_0 + s_t z_\frac{\alpha}{2}, t_0 + s_t z_{1-\frac{\alpha}{2}}\right]}
#'   This method is valid when the statistic is normally distributed around the estimate.
#' }
#' \item{`"norm"` (normal approximation)}{
#'   \deqn{\left[2t_0 - \hat{t} + s_t z_\frac{\alpha}{2}, 2t_0 - \hat{t} + s_t z_{1-\frac{\alpha}{2}}\right]}
#'
#' This involves subtracting the "bias" (\eqn{\hat{t} - t_0}) from the estimate \eqn{t_0} and using a standard Wald-type confidence interval. This method is valid when the statistic is normally distributed.
#' }
#' \item{`"basic"`}{
#'   \deqn{\left[2t_0 - t^{(1-\frac{\alpha}{2})}, 2t_0 - t^{(\frac{\alpha}{2})}\right]}
#' }
#'
#' \item{`"perc"` (percentile confidence interval)}{
#'   \deqn{\left[t^{(\frac{\alpha}{2})}, t^{(1-\frac{\alpha}{2})}\right]}
#' }
#'
#' \item{`"bc"` (bias-corrected percentile confidence interval)}{
#'   \deqn{\left[t^{(l)}, t^{(u)}\right]}
#'
#' \eqn{l = \Phi\left(2z_0 + z_\frac{\alpha}{2}\right)}, \eqn{u = \Phi\left(2z_0 + z_{1-\frac{\alpha}{2}}\right)}, where \eqn{\Phi(.)} is the normal cumulative density function (i.e., [pnorm()]) and \eqn{z_0 = \Phi^{-1}(q)} where \eqn{q} is the proportion of bootstrap estimates less than the original estimate \eqn{t_0}. This is similar to the percentile confidence interval but changes the specific quantiles of the bootstrap estimates to use, correcting for bias in the original estimate. It is described in Xu et al. (2020). When \eqn{t^0} is the median of the bootstrap distribution, the `"perc"` and `"bc"` intervals coincide.
#' }
#' \item{`"bca"` (bias-corrected and accelerated confidence interval)}{
#'   \deqn{\left[t^{(l)}, t^{(u)}\right]}
#'
#' \eqn{l = \Phi\left(z_0 + \frac{z_0 + z_\frac{\alpha}{2}}{1-a(z_0+z_\frac{\alpha}{2})}\right)}, \eqn{u = \Phi\left(z_0 + \frac{z_0 + z_{1-\frac{\alpha}{2}}}{1-a(z_0+z_{1-\frac{\alpha}{2}})}\right)}, using the same definitions as above, but with the additional acceleration parameter \eqn{a}, where \eqn{a = \frac{1}{6}\frac{\sum{L^3}}{(\sum{L^2})^{3/2}}}. \eqn{L} is the empirical influence value of each unit, which is computed using the regression method described in \pkgfun{boot}{empinf}. When \eqn{a=0}, the `"bca"` and `"bc"` intervals coincide. The acceleration parameter corrects for bias and skewness in the statistic. It can only be used when clusters are absent and the number of bootstrap replications is larger than the sample size. Note that BCa intervals cannot be requested when `simple = TRUE` and there is randomness in the `statistic` supplied to `fwb()`.
#' }
#' }
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
#' fwb_out <- fwb(infert, fit_fun, R = 199,
#'                verbose = FALSE)
#'
#' # Bias corrected percentile interval
#' bcci <- fwb.ci(fwb_out, index = "spontaneous",
#'                type = "bc")
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
  chk::chk_number(conf)
  chk::chk_range(conf, c(0, 1), inclusive = FALSE)

  chk::chk_character(type)
  type <- match_arg(type, c(.allowed_ci.types(), "all"),
                    several.ok = TRUE)

  if (any(type == "all")) {
    if (is_null(fwb.out[["cluster"]])) {
      type <- .allowed_ci.types()
    }
    else {
      type <- setdiff(.allowed_ci.types(), "bca")
    }
  }

  if (any(type == "bca")) {
    msg <- character(0L)

    if (is_not_null(fwb.out[["cluster"]])) {
      msg <- c(msg, "BCa confidence intervals cannot be used with clusters")
    }

    if (fwb.out[["R"]] <= nrow(fwb.out[["data"]])) {
      msg <- c(msg, "BCa confidence intervals cannot be computed when there are fewer bootstrap replications than units in the original dataset")
    }

    if (isTRUE(attr(fwb.out, "simple", TRUE)) &&
        isTRUE(attr(fwb.out, "random_statistic", TRUE))) {
      msg <- c(msg, 'BCa confidence intervals cannot be computed when there is randomness in `statistic` and `simple = TRUE` in the call to `fbw()`. See `vignette("fwb-rep")` for details')
    }

    if (is_not_null(msg)) {
      if (all(type == "bca")) {
        .err(msg[1L])
      }

      .wrn(msg[1L])

      type <- setdiff(type, "bca")
    }
  }

  index <- check_index(index, fwb.out[["t"]])

  t <- fwb.out[["t"]][, index, drop = FALSE]
  t0 <- fwb.out[["t0"]][index]

  if (anyNA(t)) {
    .err("some bootstrap estimates are `NA`; cannot calculate confidence intervals")
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

  output <- c(list(R = length(fins), t0 = hinv(t0), call = call),
              setNames(vector("list", length(type)), type))

  if (!identical(t, h(t))) {
    .err("`h` can only be `identity()`. Other transformations are not supported")
  }

  for (i in type) {
    output[[i]] <- compute_ci(i, t, t0, conf, 1L, hinv, fwb.out)
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

  ci.types <- intersect(names(ci.out)[lengths(ci.out) > 0L],
                        .allowed_ci.types())

  ntypes <- length(ci.types)
  t0 <- ci.out[["t0"]]

  if (is_null(hinv)) {
    use_hinv <- FALSE
    hinv <- identity
  }
  else {
    use_hinv <- TRUE
  }

  t0 <- hinv(t0)

  digs <- ceiling(log10(abs(t0)))
  digs <- {
    if (digs <= 0) 4
    else if (digs >= 4) 0
    else 4 - digs
  }

  intlabs <- setNames(character(length(ci.types)),
                      ci.types)

  basrg <- bcperg <- perg <- bcarg <- NULL

  if ("wald" %in% ci.types) {
    intlabs["wald"] <- "Wald"
  }

  if ("norm" %in% ci.types) {
    intlabs["norm"] <- "Normal"
  }

  if ("basic" %in% ci.types) {
    intlabs["basic"] <- "Basic"
    basrg <- range(ci.out[["basic"]][, 2:3])
  }

  if ("bc" %in% ci.types) {
    intlabs["bc"] <- "BC Percentile"
    bcperg <- range(ci.out[["bc"]][, 2:3])
  }

  if ("perc" %in% ci.types) {
    intlabs["perc"] <- "Percentile"
    perg <- range(ci.out[["perc"]][, 2:3])
  }

  if ("bca" %in% ci.types) {
    intlabs["bca"] <- "BCa"
    bcarg <- range(ci.out[["bca"]][, 2:3])
  }

  level <- 100 * attr(ci.out, "conf", TRUE)

  intervals <- do.call("rbind", lapply(ci.types, function(i) {
    hinv(.tail(ci.out[[i]][1L,], 2L))
  }))

  intervals_text <- apply(format(round(intervals, digs)), 1L, function(z) {
    sprintf("(%s, %s)", trimws(z[1L]), trimws(z[2L]))
  })

  text <- format(c(intlabs[ci.types], intervals_text),
                 justify = "c")

  labels <- text[seq_along(ci.types)]
  intervals_text <- text[-seq_along(ci.types)]

  level_text <- format(c("Level", sprintf("%s%%", level)),
                       justify = "c")

  R <- ci.out[["R"]]
  cat(format(c("BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS",
               sprintf("Based on %s bootstrap replicates", R)),
             justify = "c"), sep = "\n")
  cat("\nCALL :\n")
  dput(cl, control = NULL)
  cat("\nIntervals :\n")

  nc <- 2L #no. columns in output
  for (k0 in seq_len(ntypes)) {
    if ((k0 - 1) %% nc == 0) {
      upper_text <- paste(level_text[1L], labels[k0])
      lower_text <- paste(level_text[2L], intervals_text[k0])
    }
    else {
      upper_text <- paste(upper_text, labels[k0])
      lower_text <- paste(lower_text, intervals_text[k0])
    }

    if (k0 == ntypes || k0 %% nc == 0) {
      cat(sprintf("\n%s\n%s\n", upper_text, lower_text))
    }
  }

  if (is_not_null(cl[["h"]])) {
    if (is_null(cl[["hinv"]]) && !use_hinv) {
      cat("\nCalculations and Intervals on Transformed Scale\n")
    }
    else {
      cat("\nCalculations on Transformed Scale; Intervals on Original Scale\n")
    }
  }
  else if (is_null(cl[["hinv"]]) && !use_hinv) {
    cat("\nCalculations and Intervals on Original Scale\n")
  }
  else {
    cat("\nCalculations on Original Scale but Intervals Transformed\n")
  }

  if (is_not_null(basrg)) {
    if ((basrg[1L] <= 1) || (basrg[2L] >= R))
      cat("Warning : Basic Intervals used Extreme Quantiles\n")
    if ((basrg[1L] <= 10) || (basrg[2L] >= R - 9))
      cat("Some basic intervals may be unstable\n")
  }

  if (is_not_null(bcperg)) {
    if ((bcperg[1L] <= 1) || (bcperg[2L] >= R))
      cat("Warning : Bias-Corrected Percentile Intervals used Extreme Quantiles\n")
    if ((bcperg[1L] <= 10) || (bcperg[2L] >= R - 9))
      cat("Some bias-corrected percentile intervals may be unstable\n")
  }

  if (is_not_null(perg)) {
    if ((perg[1L] <= 1) || (perg[2L] >= R))
      cat("Warning : Percentile Intervals used Extreme Quantiles\n")
    if ((perg[1L] <= 10) || (perg[2L] >= R - 9))
      cat("Some percentile intervals may be unstable\n")
  }

  if (is_not_null(bcarg)) {
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

  allowed_cis <- intersect(.allowed_ci.types(),
                           names(x))

  if ("all" %in% type) {
    type <- allowed_cis
  }
  else {
    type <- match_arg(type, allowed_cis)
  }

  out <- setNames(lapply(type, function(t) {
    setNames(.tail(x[[t]][1L, ], 2L), c("L", "U"))
  }), type)

  attr(out, "conf") <- {
    if (is_not_null(attr(x, "conf", TRUE))) attr(x, "conf", TRUE)
    else x[[4L]][, 1L]
  }

  out
}

