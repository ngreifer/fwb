na.rem <- function(x) {
  #A faster na.omit for vectors
  x[!is.na(x)]
}

check_if_zero <- function(x) {
  # this is the default tolerance used in all.equal
  tolerance <- sqrt(.Machine$double.eps)
  abs(x) < tolerance
}

#Check if all values are the same
all_the_same <- function(x, na.rm = TRUE) {
  if (anyNA(x)) {
    x <- na.rem(x)
    if (!na.rm) {
      return(is_null(x))
    }
  }

  if (is.numeric(x)) check_if_zero(max(x) - min(x))
  else all(x == x[1L])
}

#Format percentage for CI labels
fmt.prc <- function(probs, digits = 3L) {
  paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), "%")
}

squish <- function(p, lo = 1e-6, hi = 1 - lo) {
  if (lo > -Inf)
    p[p < lo] <- lo

  if (hi < Inf)
    p[p > hi] <- hi

  p
}

check_index <- function(index, t, several.ok = FALSE) {
  if (is_null(index)) {
    return(1L)
  }

  if (ncol(t) == 1L) {
    if (!(isTRUE(all.equal(index, 1)) ||
          (is_not_null(colnames(t)) && isTRUE(all.equal(index, colnames(t)[1L]))))) {
      arg::wrn("only one statistic is available; ignoring {.arg index}")
    }

    return(1L)
  }

  if ((!is.character(index) && !is.numeric(index)) || is_not_null(dim(index))) {
    if (several.ok) {
      arg::err("{.arg index} must be a character or numeric vector indicating the names or indices of the desired statistics")
    }
    else {
      arg::err("{.arg index} must be a string or number indicating the name or index of the desired statistic")
    }
  }

  index <- unique(drop(index))

  if (!several.ok && length(index) > 1L) {
    arg::err("{.arg index} must have length one")
  }

  if (is.numeric(index)) {
    if (!rlang::is_integerish(index)) {
      if (several.ok) {
        arg::err("{.arg index} must be a vector of positive integers")
      }
      else {
        arg::err("{.arg index} must be a positive integer")
      }
    }

    if (any(index > ncol(t))) {
      arg::err(sprintf("{.arg index} must be between 1 and %s", ncol(t)))
    }
  }
  else {
    if (is_null(colnames(t))) {
      arg::err("the estimates don't have names, so {.arg index} must be numeric")
    }

    index <- match(index, colnames(t))

    if (anyNA(index)) {
      arg::err("all entries in {.arg index} must be the names of available statistics to compute. The following are allowed: {.val {colnames(t)}}")
    }
  }

  index
}

add_quotes <- function(x, quotes = 2L) {
  if (isFALSE(quotes)) {
    return(x)
  }

  if (isTRUE(quotes)) {
    quotes <- '"'
  }

  if (rlang::is_string(quotes)) {
    return(paste0(quotes, x, str_rev(quotes)))
  }

  if (!rlang::is_scalar_integerish(quotes) || quotes > 2 || quotes < 0) {
    stop("`quotes` must be boolean, 1, 2, or a string.")
  }

  if (quotes == 0) {
    return(x)
  }

  x <- {
    if (quotes == 1) sprintf("'%s'", x)
    else sprintf('"%s"', x)
  }

  x
}

str_rev <- function(x) {
  vapply(lapply(strsplit(x, NULL), rev), paste, character(1L), collapse = "")
}

.tail <- function(x, n = 1L) {
  arg::arg_count(n)
  arg::arg_gt(n, 0)

  l <- length(x)
  x[seq(max(1L, l - n + 1L), l)]
}

.attr <- function(x, which, exact = TRUE) {
  attr(x, which, exact = exact)
}

is_null <- function(x) {length(x) == 0L}
is_not_null <- function(x) {!is_null(x)}
`%or%` <- function(x, y) {
  # like `%||%` but works for non-NULL length 0 objects
  if (is_null(x)) y else x
}

.set_class <- function(x, .class = NULL) {
  if (is_null(.class)) {
    return(x)
  }

  class(x) <- c(class(x), .class)

  x
}

with_seed_preserved <- function(expr, new_seed = NULL) {

  old_seed <- list(random_seed = get(".Random.seed", globalenv(), mode = "integer",
                                     inherits = FALSE),
                   rng_kind = RNGkind())

  if (is_null(old_seed)) {
    on.exit({
      set.seed(seed = NULL)
      rm(".Random.seed", envir = globalenv())
    }, add = TRUE)
  }
  else {
    on.exit({
      .RNGkind <- get("RNGkind")
      .RNGkind(old_seed$rng_kind[[1L]], normal.kind = old_seed$rng_kind[[2L]])
      sample_kind <- old_seed$rng_kind[[3L]]

      if (identical(sample_kind, "Rounding")) {
        suppressWarnings(.RNGkind(sample.kind = sample_kind))
      }
      else {
        .RNGkind(sample.kind = sample_kind)
      }

      assign(".Random.seed", old_seed$random_seed, globalenv())
    }, add = TRUE)
  }

  if (is_not_null(new_seed)) {
    assign(".Random.seed", value = new_seed, envir = globalenv())
  }

  expr
}

guess_num_workers <- function(cl = NULL) {
  if (is_null(cl)) {
    return(1L)
  }

  if (identical(cl, "future")) {
    rlang::check_installed("future")
    return(future::nbrOfWorkers())
  }

  if (inherits(cl, "cluster")) {
    return(length(cl))
  }

  if (is.numeric(cl)) {
    if (.Platform$OS.type == "windows") {
      return(1L)
    }

    return(max(1L, as.integer(cl)))
  }

  1L
}
