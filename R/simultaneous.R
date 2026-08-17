# Simultaneous CIs and p-values

# For each bootstrap draw, the narrowest pointwise level at which it enters the band.
#
# The pointwise interval for parameter `j` at level `l` is
# `[norm_inter(t_j, p), norm_inter(t_j, 1 - p)]` with `p = (1 - l) / 2`. `norm_inter()` is
# increasing in `p` and returns the `i`th order statistic exactly at `p = i / (R + 1)`, so
# writing `i` for the rank of `t[r, j]` within column `j`:
#
#     t[r, j] >= lower  <=>  p <= i / (R + 1)
#     t[r, j] <= upper  <=>  p <= 1 - i / (R + 1)
#
# The coordinate is therefore inside its own interval exactly when
# `p <= min(i, R + 1 - i) / (R + 1)`, and the draw is inside the whole band when that holds
# for every `j` -- i.e. at the smallest such `p` across parameters.
#
# Ties are broken so that whichever endpoint binds is the one used: the highest tied rank
# against the lower endpoint, the lowest against the upper.
sup_t_entry_levels <- function(t) {
  R <- nrow(t)

  lo <- apply(t, 2L, rank, ties.method = "min")
  hi <- apply(t, 2L, rank, ties.method = "max")

  p <- pmin(hi, R + 1L - lo) / (R + 1)

  1 - 2 * apply(p, 1L, min)
}

# The sup-t percentile level: the narrowest pointwise level whose intervals jointly contain
# at least `level` of the bootstrap draws.
#
# Coverage at level `l` is the proportion of draws whose entry level is at most `l`, so this
# needs no search -- it is the `ceiling(level * R)`th smallest entry level. Two ranks per
# parameter and one partial sort, with no interval ever constructed.
#
# The result is nudged a few ulps wider, which is not cosmetic. At the exact level the
# binding draw sits *on* an endpoint, and `compute_ci()` does not receive `p` -- it receives
# the level and recomputes `p = (1 - level) / 2`. Since the level is near 1, that round trip
# loses about an ulp of `p`, enough to flip `floor((R + 1) * p)` from `i` to `i - 1`; the
# endpoint then lands a hair inside the binding draw and coverage comes out one draw short.
# The nudge puts `p` safely on the correct side of the floor, and moves each endpoint by
# ~1e-16 of the gap between adjacent order statistics -- unchanged to full double precision.
sup_t_level <- function(t, level) {
  entry <- sup_t_entry_levels(t)

  m <- ceiling(level * nrow(t))

  sort(entry, partial = m)[m] + 16 * .Machine$double.eps
}

simultaneous_ci_level <- function(object, level = .95, index = seq_len(ncol(object[["t"]])), ci.type = "perc") {

  arg::arg_supplied(object)
  arg::arg_is(object, "fwb")

  arg::arg_number(level)
  arg::arg_between(level, c(0, 1), inclusive = FALSE)

  index <- check_index(index, object[["t"]], several.ok = TRUE)

  k <- length(index)

  if (k == 1L) {
    return(level)
  }

  arg::arg_string(ci.type)

  if (!ci.type %in% c("perc", "wald")) {
    arg::err("simultaneous inference can only be used with Wald or percentile intervals")
  }

  if (ci.type == "perc") {
    t <- object[["t"]][, index, drop = FALSE]

    if (!all(is.finite(t))) {
      arg::err("simultaneous inference cannot be used when some bootstrap estimates are {.val {NA}} or non-finite")
    }

    new_level <- sup_t_level(t, level)
  }
  else {
    rlang::check_installed("mvtnorm")

    arg::arg_between(level, c(0, 1), inclusive = FALSE)

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

simultaneous_p_value <- function(object, p.values, index = seq_len(ncol(object[["t"]])), ci.type = "perc") {
  arg::arg_supplied(object)
  arg::arg_is(object, "fwb")

  arg::arg_supplied(p.values)
  arg::arg_numeric(p.values)
  arg::arg_between(p.values, c(0, 1))

  index <- check_index(index, object[["t"]], several.ok = TRUE)

  arg::arg_string(ci.type)

  k <- length(index)

  if (k == 1L) {
    return(p.values)
  }

  if (!ci.type %in% c("perc", "wald")) {
    arg::err("simultaneous inference can only be used with Wald or percentile intervals")
  }

  if (ci.type == "perc") {
    t <- object[["t"]][, index, drop = FALSE]

    if (!all(is.finite(t))) {
      arg::err("simultaneous inference cannot be used when some bootstrap estimates are {.val {NA}} or non-finite")
    }

    #The simultaneous p-value is one minus the coverage of the band whose pointwise level
    #is `1 - p`. Once the entry levels are in hand that is a comparison, not an interval
    #calculation -- and it is the same quantity `simultaneous_ci_level()` inverts, computed
    #the same way, which is what keeps the p-values consistent with the bands.
    entry <- sup_t_entry_levels(t)

    p <- vapply(1 - p.values, function(level) mean(entry > level), numeric(1L))
  }
  else {
    rlang::check_installed("mvtnorm")

    v <- cov(object[["t"]][, index, drop = FALSE])

    zeros <- check_if_zero(diag(v))

    v <- cov2cor(v[!zeros, !zeros, drop = FALSE])

    z <- abs(qnorm(p.values / 2))

    p <- rep.int(0.0, length(index))

    p[!zeros] <- vapply(z[!zeros], function(zi) {
      1 - mvtnorm::pmvnorm(lower = rep.int(-zi, sum(!zeros)),
                           upper = rep.int(zi, sum(!zeros)),
                           corr = v,
                           abseps = 1e-5,
                           maxpts = 1e6)
    }, numeric(1L))
  }

  p
}
