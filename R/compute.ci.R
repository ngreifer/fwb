#Functions for computing and inverting confidence intervals

# The fewest bootstrap replications each interval type can be computed from
min_R_for_ci <- function(type) {
  if (type == "cheap") {
    return(1L)
  }

  2L
}

# Whether `type` can be computed at all from these replicates, with a message naming the
# type and its requirement when it cannot.
check_ci_feasible <- function(type, R, boot.out = NULL) {
  min_R <- min_R_for_ci(type)

  if (R < min_R) {
    arg::err("the {.val {type}} confidence interval requires at least {min_R} bootstrap replications; there {?is/are} only {R}")
  }

  if (identical(type, "bca") && is_not_null(boot.out)) {
    n <- nrow(boot.out[["data"]])

    if (R <= n) {
      arg::err("the {.val {type}} confidence interval requires more bootstrap replications than units in the original dataset; there are {R} replications and {n} units")
    }

    if (is_not_null(boot.out[["cluster"]])) {
      arg::err("the {.val {type}} confidence interval cannot be used with clusters")
    }
  }

  invisible(NULL)
}

# Compute confidence intervals
compute_ci <- function(type, t, t0, conf = .95, index = 1, hinv = identity, boot.out = NULL) {
  check_ci_feasible(type, nrow(t), boot.out)

  alpha <- (1 + c(-conf, conf)) / 2

  if (type %in% c("wald", "norm")) {
    bias <- switch(type,
                   wald = 0,
                   norm = colMeans(t[, index, drop = FALSE]) - t0[index])

    sds <- apply(t[, index, drop = FALSE], 2L, sd)

    zcrit <- abs(qnorm((1 - conf) / 2))

    out <- cbind(conf,
                 hinv(t0[index] - bias - sds * zcrit),
                 hinv(t0[index] - bias + sds * zcrit))

    return(out)
  }

  if (type == "cheap") {
    sds <- sqrt(colMeans(sweep(t[, index, drop = FALSE], 2L, t0[index])^2))

    tcrit <- abs(qt((1 - conf) / 2, df = nrow(t)))

    out <- cbind(conf,
                 hinv(t0[index] - sds * tcrit),
                 hinv(t0[index] + sds * tcrit))

    return(out)
  }

  if (type == "basic") {
    out <- do.call("rbind", lapply(index, function(i) {
      qq <- norm_inter(t[, i], rev(alpha))

      cbind(conf,
            matrix(qq[, 1L], ncol = 2L),
            matrix(hinv(2 * t0[i] - qq[, 2L]), ncol = 2L))
    }))

    return(out)
  }

  if (type == "perc") {
    adj.alpha <- function(i) {
      alpha
    }
  }
  else {
    s <- vapply(index, function(i) {
      min(max(1, sum(t[, i] < t0[i])), nrow(t) - 1)
    }, numeric(1L))

    w <- qnorm(s / nrow(t))

    zalpha <- qnorm(alpha)

    if (type == "bca") {
      if (is_null(boot.out)) {
        arg::err("BCa confidence intervals cannot be computed")
      }

      if (is_not_null(boot.out[["cluster"]])) {
        arg::err("the BCa confidence interval cannot be used with clusters")
      }

      L <- empinf.reg(boot.out, t = t[, index, drop = FALSE])

      a <- colSums(L^3, na.rm = TRUE) / (6 * colSums(L^2, na.rm = TRUE)^1.5)

      if (!all(is.finite(a))) {
        arg::err("estimated adjustment {.var a} is {.val {NA}}")
      }
    }
    else { # type == "bc"
      a <- rep.int(0, length(index))
    }

    adj.alpha <- function(i) {
      pnorm(w[i] + (w[i] + zalpha) / (1 - a[i] * (w[i] + zalpha)))
    }
  }

  do.call("rbind", lapply(seq_along(index), function(i) {
    qq <- norm_inter(t[, index[i]], adj.alpha(i))

    cbind(conf,
          matrix(qq[, 1L], ncol = 2L),
          matrix(hinv(qq[, 2L]), ncol = 2L))
  }))
}

# Normal interpolation; similar to quantile(ti, probs = alpha)
norm_inter <- function(ti, alpha = c(.025, .975)) {
  arg::arg_numeric(alpha)

  ti <- ti[is.finite(ti)]
  R <- length(ti)

  # Scaled quantile ranks
  scaled_idx <- (R + 1) * alpha
  k <- floor(scaled_idx)

  if (any(k < 1) || any(k >= R)) {
    arg::wrn("extreme order statistics used as endpoints")

    k[k < 1] <- 1
    k[k >= R] <- R - 1
  }

  # Corresponding z-scores
  z_k     <- qnorm(k / (R + 1))
  z_k1    <- qnorm((k + 1) / (R + 1))
  z_alpha <- qnorm(squish(alpha, 1e-9))

  # Interpolation weight on z-scale
  w <- squish((z_alpha - z_k) / (z_k1 - z_k), 0, 1)

  if (is.unsorted(ti)) {
    ti <- sort(ti, partial = sort(unique(c(k, k + 1L))))
  }

  # Final interpolated quantile
  q_interp <- (1 - w) * ti[k] + w * ti[k + 1L]

  cbind(round(scaled_idx, 2L), q_interp)
}

# Regression-based empirical influence function. Centered coefficients from
# regression of estimates on bootstrap weights. Yields 1 value per unit per
# quantity of interest. Slow due to R x n design matrix -- often slower than the
# `fwb()` run that produced `t`.
#
# `crossprod()`-based normal equations are the obvious replacement and were
# measured: 30x faster under an optimized BLAS, but no faster (and up to 1.7x
# slower when R is close to n) under the reference BLAS that R ships by default,
# since `.lm.fit()` uses level-1 LINPACK and so does not care which BLAS is
# linked. They also lose 5 digits at `R = n + 1`, which `fwb.ci()` permits.
# Not worth it; see dev/review-2026-08.md.
empinf.reg <- function(boot.out, t = NULL) {
  n <- NROW(boot.out[["data"]])
  f <- fwb.array(boot.out)

  if (is_null(t)) {
    t <- boot.out[["t"]]
  }

  if (all(is.finite(t))) {
    out <- as.matrix(.lm.fit(x = f / n, y = t)$coefficients)
    out[] <- sweep(out, 2L, colMeans(out))
  }
  else {
    out <- matrix(NA_real_, nrow = n, ncol = ncol(t))

    for (i in seq_len(ncol(t))) {
      fins <- which(is.finite(t[, i]))
      ti <- t[fins, i]

      l <- .lm.fit(x = f[fins, , drop = FALSE] / n, y = ti)$coefficients
      out[fins, i] <- l - mean(l)
    }
  }

  colnames(out) <- colnames(t)

  out
}

#Invert confidence intervals to get p-values
invert_ci <- function(type, t, t0, null = 0, index = 1L, h = identity, boot.out = NULL) {
  #A p-value here is the level at which the interval's endpoint lands on `null`, so it
  #needs exactly what the interval needs.
  check_ci_feasible(type, nrow(t), boot.out)

  if (type %in% c("wald", "norm")) {
    bias <- switch(type,
                   wald = 0,
                   norm = colMeans(t[, index, drop = FALSE]) - t0[index])

    sds <- apply(t[, index, drop = FALSE], 2L, sd)

    alpha <- pnorm(-abs(t0[index] - bias - h(null)) / sds)

    p <- 2 * alpha

    return(p)
  }

  if (type == "cheap") {
    sds <- sqrt(colMeans(sweep(t[, index, drop = FALSE], 2L, t0[index])^2))

    alpha <- pt(-abs(t0[index] - h(null)) / sds, df = nrow(t))

    p <- 2 * alpha

    return(p)
  }

  if (type == "basic") {
    alpha <- vapply(index, function(i) {
      norm_inter_invert(t[, i], 2 * t0[i] - h(null))
    }, numeric(1L))

    p <- 2 * pmin(alpha, 1 - alpha)

    return(p)
  }

  quant <- apply(t[, index, drop = FALSE], 2L, norm_inter_invert, h(null))

  if (type == "perc") {
    alpha <- quant
  }
  else if (type == "bc") {
    s <- vapply(index, function(i) {
      min(max(1, sum(t[, i] < t0[i])), nrow(t) - 1)
    }, numeric(1L))

    w <- qnorm(s / nrow(t))

    alpha <- pnorm(qnorm(quant) - 2 * w)
  }
  else if (type == "bca") {
    if (is_null(boot.out)) {
      arg::err("BCa p-values cannot be computed")
    }

    if (is_not_null(boot.out[["cluster"]])) {
      arg::err("BCa p-values cannot be used with clusters")
    }

    L <- empinf.reg(boot.out, t = t[, index, drop = FALSE])

    a <- colSums(L^3, na.rm = TRUE) / (6 * colSums(L^2, na.rm = TRUE)^1.5)

    if (!all(is.finite(a))) {
      arg::err("estimated adjustment {.var a} is {.val {NA}}")
    }

    s <- vapply(index, function(i) {
      min(max(1, sum(t[, i] < t0[i])), nrow(t) - 1)
    }, numeric(1L))

    w <- qnorm(s / nrow(t))

    z <- qnorm(quant) - w

    alpha <- rep.int(0.0, length(index))

    fins <- is.finite(z)

    if (any(fins)) {
      alpha[fins] <- pnorm(z[fins] / (1 + z[fins] * a[fins]) - w[fins])
    }
  }

  2 * pmin(alpha, 1 - alpha)
}

# Normal interpolation inverse; similar to mean(ti <= null)
norm_inter_invert <- function(ti, null = 0) {
  arg::arg_number(null)

  ti <- ti[is.finite(ti)]
  R <- length(ti)

  if (is.unsorted(ti)) {
    ti <- sort(ti)
  }

  if (null < ti[1L]) {
    return(0)
  }

  if (null >= ti[R]) {
    return(1)
  }

  if (getRversion() >= "4.5.0") {
    k <- findInterval(null, ti, checkSorted = FALSE, checkNA = FALSE)
  }
  else {
    k <- findInterval(null, ti)
  }

  # Exact quantile if null equal to observed value
  if (ti[k] == null) {
    return(unname(k / (R + 1L)))
  }

  r <- (null - ti[k]) / (ti[k + 1L] - ti[k])

  # Normal quantiles for interpolation
  z_k   <- qnorm(k / (R + 1L))
  z_k1  <- qnorm((k + 1L) / (R + 1L))

  # Interpolate on z-scale
  z_interp <- z_k + r * (z_k1 - z_k)

  unname(pnorm(z_interp))
}

.allowed_ci.types <- function() {
  c("perc", "bc", "wald", "norm", "cheap", "basic", "bca")
}
