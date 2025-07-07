# Simultaneous CIs and p-values
simultaneous_ci_level <- function(object, level = .95, index = 1:ncol(object[["t"]]), ci.type = "perc") {

  chk::chk_not_missing(object, "`object`")
  chk::chk_is(object, "fwb")

  chk::chk_number(level)
  chk::chk_range(level, c(0, 1), inclusive = FALSE)

  index <- check_index(index, object[["t"]], several.ok = TRUE)

  k <- length(index)

  if (k == 1L) {
    return(level)
  }

  chk::chk_string(ci.type)

  if (!ci.type %in% c("perc", "wald")) {
    .err("simultaneous inference can only be used with Wald or percentile intervals")
  }

  if (ci.type == "perc") {
    estimates <- t(object[["t"]][, index, drop = FALSE])
    R <- ncol(estimates)

    fun <- function(l) {
      suppressWarnings({
        ci.out <- compute_ci(ci.type, t = object[["t"]], t0 = object[["t0"]],
                             conf = l, index = index, boot.out = object)
      })

      nc <- ncol(ci.out)

      interval <- ci.out[, c(nc - 1L, nc), drop = FALSE]

      all_est_above <- colSums(estimates >= interval[, 1L]) == k
      all_est_below <- colSums(estimates <= interval[, 2L]) == k

      coverage <- mean(all_est_above & all_est_below)

      R * abs(coverage - level) + l
    }

    new_level <- optimize(fun, interval = c(level, 1 - (1 - level) / k),
                          tol = 1e-10)$minimum
  }
  else {
    rlang::check_installed("mvtnorm")

    chk::chk_range(level, c(.5, 1), inclusive = FALSE)

    v <- cov(object[["t"]][, index, drop = FALSE])

    zeros <- check_if_zero(diag(v))

    v <- cov2cor(v[!zeros, !zeros, drop = FALSE])

    z_crit <- mvtnorm::qmvnorm(level,
                               tail = "both.tails",
                               corr = v, keepAttr = FALSE,
                               ptol = .0001, maxiter = 1e4)$quantile

    new_level <- 1 - 2 * pnorm(-abs(z_crit))
  }

  new_level
}

simultaneous_p_value <- function(object, p.values, index = 1:ncol(object[["t"]]), ci.type = "perc") {
  chk::chk_not_missing(object, "`object`")
  chk::chk_is(object, "fwb")

  chk::chk_not_missing(p.values, "`p.values`")
  chk::chk_numeric(p.values)
  chk::chk_range(p.values, c(0, 1))

  index <- check_index(index, object[["t"]], several.ok = TRUE)

  chk::chk_string(ci.type)

  k <- length(index)

  if (k == 1L) {
    return(p.values)
  }

  if (!ci.type %in% c("perc", "wald")) {
    .err("simultaneous inference can only be used with Wald or percentile intervals")
  }

  if (ci.type == "perc") {
    estimates <- t(object[["t"]][, index, drop = FALSE])

    p <- vapply(1 - p.values, function(level) {
      suppressWarnings({
        ci.out <- compute_ci(ci.type, t = object[["t"]], t0 = object[["t0"]],
                             conf = level, index = index, boot.out = object)
      })

      nc <- ncol(ci.out)

      interval <- ci.out[, c(nc - 1L, nc), drop = FALSE]

      all_est_above <- colSums(estimates >= interval[, 1L]) == k
      all_est_below <- colSums(estimates <= interval[, 2L]) == k

      coverage <- mean(all_est_above & all_est_below)

      1 - coverage
    }, numeric(1L))
  }
  else {
    rlang::check_installed("mvtnorm")

    v <- cov(object[["t"]][, index, drop = FALSE])

    zeros <- check_if_zero(diag(v))

    v <- cov2cor(v[!zeros, !zeros, drop = FALSE])

    z <- abs(qnorm(p.values / 2))

    p <- rep.int(0.0, length(index))

    p[!zeros] <- vapply(z, function(zi) {
      1 - mvtnorm::pmvnorm(lower = rep.int(-zi, sum(!zeros)),
                           upper = rep.int(zi, sum(!zeros)),
                           corr = v,
                           abseps = 1e-5,
                           maxpts = 1e6)
    }, numeric(1L))
  }

  p
}

