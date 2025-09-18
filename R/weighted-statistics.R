#' Calculate weighted statistics
#'
#' These functions are designed to compute weighted statistics (mean, variance, standard deviation, covariance, correlation, median and quantiles) to perform weighted transformation (scaling, centering, and standardization) using bootstrap weights. These automatically extract bootstrap weights when called from within [fwb()] to facilitate computation of bootstrap statistics without needing to think to hard about how to correctly incorporate weights.
#'
#' @param x a numeric variable; for `w_cov()` and `w_cor()`, a numeric matrix.
#' @param w optional; a numeric vector of weights with length equal to the length of `x` (or number of rows if `x` is a matrix). If unspecified, will use bootstrapped weights when called from within [fwb()] or [vcovFWB()] and unit weights (i.e., for unweighted estimates) otherwise.
#' @param na.rm logical; whether to exclude `NA` values in the weights and `x` when computing statistics. Default is `FALSE` for the weighted statistics (like with their unweighted counterparts) and `TRUE` for weighted transformations.
#' @param probs,names,type,digits see [quantile()]. Only `type = 7` is allowed.
#' @param scale,center logical; whether to scale or center the variable.
#'
#' @return `w_mean()`, `w_var()`, `w_sd()`, and `w_median()` return a numeric scalar. `w_cov()` and `w_cor()` return numeric matrices. `w_quantile()` returns a numeric vector of length equal to `probs`. `w_std()`, `w_scale()`, and `w_center()` return numeric vectors of length equal to the length of `x`.
#'
#' @details
#' These function automatically incorporate bootstrap weights when used inside [fwb()] or [vcovFWB()]. This works because `fwb()` and `vcovFWB()` temporarily set an option with the environment of the function that calls the estimation function with the sampled weights, and the `w_*()` functions access the bootstrap weights from that environment, if any. So, using, e.g., `w_mean()` inside the function supplied to the `statistic` argument of `fwb()`, computes the weighted mean using the bootstrap weights. Using these functions outside `fwb()` works just like other functions that compute weighted statistics: when `w` is supplied, the statistics are weighted, and otherwise they are unweighted.
#'
#' See below for how each statistic is computed.
#'
#' ## Weighted statistics
#' For all weighted statistics, the weights are first rescaled to sum to 1. \eqn{w} in the formulas below refers to these weights.
#'
#' \describe{
#'   \item{`w_mean()`}{Calculates the weighted mean as
#'                     \deqn{\bar{x}_w = \sum_i{w_i x_i}}
#'                     This is the same as [weighted.mean()].}
#'   \item{`w_var()`}{Calculates the weighted variance as
#'                    \deqn{s^2_{x,w} = \frac{\sum_i{w_i(x_i - \bar{x}_w)^2}}{1 - \sum_i{w_i^2}}}}
#'   \item{`w_sd()`}{Calculates the weighted standard deviation as
#'                  \deqn{s_{x,w} = \sqrt{s^2_{x,w}}}}
#'   \item{`w_cov()`}{Calculates the weighted covariances as
#'                    \deqn{s_{x,y,w} = \frac{\sum_i{w_i(x_i - \bar{x}_w)(y_i - \bar{y}_w)}}{1 - \sum_i{w_i^2}}}
#'                    This is the same as [cov.wt()].}
#'   \item{`w_cor()`}{Calculates the weighted correlation as
#'                    \deqn{r_{x,y,w} = \frac{s_{x,y,w}}{s_{x,w}s_{y,w}}}
#'                    This is the same as [cov.wt()].}
#'   \item{`w_quantile()`}{Calculates the weighted quantiles using linear interpolation of the weighted cumulative density function.}
#'   \item{`w_median()`}{Calculates the weighted median as `w_quantile(., probs = .5)`.}
#' }
#'
#' ## Weighted transformations
#'
#' Weighted transformations use the weighted mean and/or standard deviation to re-center or re-scale the given variable. In the formulas below, \eqn{x} is the original variable and \eqn{w} are the weights.
#'
#' \describe{
#'   \item{`w_center()`}{Centers the variable at its (weighted) mean.
#'                       \deqn{x_{c,w} = x - \bar{x}_w}}
#'   \item{`w_scale()`}{Scales the variable by its (weighted) standard deviation.
#'                      \deqn{x_{s,w} = x / s_{x,w}}}
#'   \item{`w_std()`}{Centers and scales the variable by its (weighted) mean and standard deviation.
#'                    \deqn{x_{z,w} = (x - \bar{x}_w) / s_{x,w}}}
#' }
#'
#' `w_scale()` and `w_center()` are efficient wrappers to `w_std()` with `center = FALSE` and `scale = FALSE`, respectively.
#'
#' @seealso
#' * [mean()] and [weighted.mean()] for computing the unweighted and weighted mean
#' * [var()] and [sd()] for computing the unweighted variance and standard deviation
#' * [median()] and [quantile()] for computing the unweighted median and quantiles
#' * [cov()], [cor()], and [cov.wt()] for unweighted and weighted covariance and correlation matrices
#' * [scale()] for standardizing variables using arbitrary (by default, unweighted) centering and scaling factors
#'
#' @examplesIf rlang::is_installed(c("cobalt", "lmtest"))
#' # G-computation of average treatment effects using lalonde
#' # dataset
#'
#' data("lalonde", package = "cobalt")
#'
#' ate_est <- function(data, w) {
#'   fit <- lm(re78 ~ treat * (age + educ + married + race +
#'                               nodegree + re74 + re75),
#'             data = data, weights = w)
#'
#'   p0 <- predict(fit, newdata = transform(data, treat = 0))
#'   p1 <- predict(fit, newdata = transform(data, treat = 1))
#'
#'   # Weighted means using bootstrap weights
#'   m0 <- w_mean(p0)
#'   m1 <- w_mean(p1)
#'
#'   c(m0 = m0, m1 = m1, ATE = m1 - m0)
#' }
#'
#' set.seed(123, "L'Ecuyer-CMRG")
#' boot_est <- fwb(lalonde, statistic = ate_est,
#'                 R = 199, verbose = FALSE)
#' summary(boot_est)
#'
#' # Using `w_*()` data transformations inside a model
#' # supplied to vcovFWB():
#' fit <- lm(re78 ~ treat * w_center(age),
#'           data = lalonde)
#'
#' lmtest::coeftest(fit, vcov = vcovFWB, R = 500)

#' @export
w_mean <- function(x, w = NULL, na.rm = FALSE) {
  if (is_null(w)) {
    w <- .get_internal_w()
  }

  chk::chk_flag(na.rm)

  if (is_null(w)) {
    return(mean(x, na.rm = na.rm))
  }

  if (na.rm && (anyNA(x) || anyNA(w))) {
    i <- !is.na(x) & !is.na(w)
    x <- x[i]
    w <- w[i]
  }

  w <- w / sum(w)

  sum(x * w)
}

#' @export
#' @rdname w_mean
w_var <- function(x, w = NULL, na.rm = FALSE) {
  if (is_null(w)) {
    w <- .get_internal_w()
  }

  chk::chk_flag(na.rm)

  if (is_null(w)) {
    return(var(x, na.rm = na.rm))
  }

  if (na.rm && (anyNA(x) || anyNA(w))) {
    i <- !is.na(x) & !is.na(w)
    x <- x[i]
    w <- w[i]
  }

  w <- w / sum(w)

  mu <- sum(x * w)

  sum(w * (x - mu)^2) / (1 - sum(w^2))
}

#' @export
#' @rdname w_mean
w_sd <- function(x, w = NULL, na.rm = FALSE) {
  sqrt(w_var(x = x, w = w, na.rm = na.rm))
}

#' @export
#' @rdname w_mean
w_cov <- function(x, w = NULL, na.rm = FALSE) {
  if (is_null(w)) {
    w <- .get_internal_w()
  }

  x <- as.matrix(x)

  if (is_null(w)) {
    return(cov(x, use = if (na.rm) "pairwise.complete.obs" else "everything"))
  }

  if (na.rm && anyNA(w)) {
    w[is.na(w)] <- 0
  }

  if (!na.rm || !anyNA(x)) {
    w <- w / sum(w)

    x_cen <- sweep(x, 2L, colSums(w * x), check.margin = FALSE)

    return(crossprod(x_cen * w, x_cen) / (1 - sum(w^2)))
  }

  m <- matrix(NA_real_, nrow = ncol(x), ncol = ncol(x),
              dimnames = list(colnames(x), colnames(x)))

  for (i in seq_len(ncol(x))) {
    ki <- !is.na(x[, i])

    for (j in seq_len(i - 1L)) {
      kij <- which(ki & !is.na(x[, j]))

      wk <- w[kij] / sum(w[kij])

      mu_i <- sum(wk * x[kij, i])
      mu_j <- sum(wk * x[kij, j])

      m[i, j] <- sum(wk * (x[kij, i] - mu_i) * (x[kij, j] - mu_j)) / (1 - sum(wk^2))
      m[j, i] <- m[i, j]
    }

    ki <- which(ki)
    wk <- w[ki] / sum(w[ki])

    mu_i <- sum(wk * x[ki, i])

    m[i, i] <- sum(wk * (x[ki, i] - mu_i)^2) / (1 - sum(wk^2))
  }

  m
}

#' @export
#' @rdname w_mean
w_cor <- function(x, w = NULL) {
  if (is_null(w)) {
    w <- .get_internal_w()
  }

  x <- as.matrix(x)

  if (is_null(w)) {
    return(cor(x))
  }

  cov2cor(w_cov(x, w))
}

#' @export
#' @rdname w_mean
w_quantile <- function(x, w = NULL, probs = seq(0, 1, by = 0.25), na.rm = FALSE, names = TRUE,
                       type = 7L, digits = 7L) {

  if (is_null(w)) {
    w <- .get_internal_w()
  }

  chk::chk_numeric(probs)
  chk::chk_range(probs)

  chk::chk_flag(na.rm)
  chk::chk_flag(names)

  chk::chk_equal(type, 7L)

  if (names) {
    chk::chk_count(digits)
  }

  if (is_null(w) || sum(w > 0, na.rm = TRUE) <= 1) {
    return(quantile(x, probs = probs, na.rm = na.rm, names = names,
                    type = type, digits = digits))
  }

  if (na.rm && (anyNA(x) || anyNA(w))) {
    i <- !is.na(x) & !is.na(w)
    x <- x[i]
    w <- w[i]
  }

  w_0 <- w == 0

  if (any(w_0)) {
    w <- w[!w_0]
    x <- x[!w_0]
  }

  w <- w / sum(w)

  if (is.unsorted(x)) {
    x_order <- order(x)
    x <- x[x_order]
    w <- w[x_order]
  }

  p_k <- (cumsum(w) - w) / (1 - w)

  qu <- approxfun(p_k, x, rule = 2L, ties = "ordered")(probs)

  if (names && is_not_null(qu)) {
    names(qu) <- paste0(formatC(probs * 100, format = "fg",
                                width = 1L, digits = digits), "%")
    names(qu)[is.na(probs)] <- ""
  }

  qu
}

#' @export
#' @rdname w_mean
w_median <- function(x, w = NULL, na.rm = FALSE) {
  if (is_null(w)) {
    w <- .get_internal_w()
  }

  chk::chk_flag(na.rm)

  if (is_null(w)) {
    return(median(x, na.rm = na.rm))
  }

  w_quantile(x, w = w, probs = .5, na.rm = na.rm, names = FALSE)
}

#' @export
#' @rdname w_mean
w_std <- function(x, w = NULL, na.rm = TRUE, scale = TRUE, center = TRUE) {
  if (is_null(w)) {
    w <- .get_internal_w()
  }

  chk::chk_flag(scale)
  chk::chk_flag(center)

  if (is_null(w)) {
    mu <- {
      if (center) mean(x, na.rm = na.rm)
      else 0
    }

    v <- {
      if (scale) var(x, na.rm = na.rm)
      else 1
    }

    return((x - mu) / sqrt(v))
  }

  if (na.rm && (anyNA(x) || anyNA(w))) {
    i <- !is.na(x) & !is.na(w)
    x <- x[i]
    w <- w[i]
  }

  w <- w / sum(w)

  w_0 <- w == 0

  if (any(w_0)) {
    w <- w[!w_0]
    x <- x[!w_0]
  }

  mu <- {
    if (center) sum(w * x)
    else 0
  }

  v <- {
    if (scale) sum(w * (x - mu)^2) / (1 - sum(w^2))
    else 1
  }

  (x - mu) / sqrt(v)
}

#' @export
#' @rdname w_mean
w_scale <- function(x, w = NULL, na.rm = TRUE) {
  w_std(x, w = w, na.rm = na.rm,
        scale = TRUE, center = FALSE)
}

#' @export
#' @rdname w_mean
w_center <- function(x, w = NULL, na.rm = TRUE) {
  w_std(x, w = w, na.rm = na.rm,
        scale = FALSE, center = TRUE)
}

.get_internal_w <- function() {
  env <- getOption("fwb_internal_w_env")

  if (is_null(env)) {
    return(NULL)
  }

  .wi <- env[[".wi"]]

  if (is_null(.wi) || !inherits(.wi, "fwb_internal_w")) {
    return(NULL)
  }

  .wi
}
