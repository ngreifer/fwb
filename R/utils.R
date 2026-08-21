na.rem <- function(x) {
  #A faster na.omit for vectors
  x[!is.na(x)]
}

check_if_zero <- function(x, tolerance = sqrt(.Machine$double.eps)) {
  # this is the default tolerance used in all.equal
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

#Row-wise tabulation: `out[r, u]` is the number of times `u` appears in row `r`
#of the index matrix `i`.
#
#Two strategies, because their costs cross over. Adding `(u - 1) * R` to each draw
#puts the count for `out[r, u]` at position `r + (u - 1) * R`, which is exactly
#where a column-major `R x n` matrix keeps it -- so the whole matrix can be counted
#by one `tabulate()` call, with no per-row overhead and no transpose. That wins by a
#wide margin while `n` is small (~7x at `n = 20`, ~2.7x at `n = 64`), because there
#the per-row overhead is most of the work. But it scatters writes across `R * n`
#bins, and past `n` of roughly 400 that costs more than the overhead it saves.
#Above the threshold, count each replicate separately, transposing first so that
#each replicate's draws are contiguous in memory.
#
#The threshold predicts on `n` alone, not on `R * n`: the first branch still wins
#at `n = 100, R = 50000` (five million bins), so this is about per-row overhead
#rather than about the bin array fitting in cache.
#
#Either way the counts depend only on the multiset of each row of `i`, never on the
#order within a row, so `wtype = "multinom"` stays identical to
#`boot:::ordinary.array()`.
tabulate_rows <- function(i, n) {
  R <- nrow(i)

  if (n <= 256L && R <= .Machine$integer.max %/% n) {
    out <- tabulate(rep.int(seq_len(R), n) + (i - 1L) * R, R * n)

    dim(out) <- c(R, n)

    return(out)
  }

  it <- t(i)

  t(vapply(seq_len(R), function(r) tabulate(it[, r], n), integer(n)))
}

check_index <- function(index, t, several.ok = FALSE, .arg_index = "index") {
  if (is_null(index)) {
    return(1L)
  }

  if (ncol(t) == 1L) {
    if (!(isTRUE(all.equal(index, 1)) ||
          (is_not_null(colnames(t)) && isTRUE(all.equal(index, colnames(t)[1L]))))) {
      arg::wrn("only one statistic is available; ignoring {.arg {(.arg_index)}}")
    }

    return(1L)
  }

  if ((!is.character(index) && !is.numeric(index)) || is_not_null(dim(index))) {
    if (several.ok) {
      arg::err("{.arg {(.arg_index)}} must be a character or numeric vector indicating the names or indices of the desired statistics")
    }
    else {
      arg::err("{.arg {(.arg_index)}} must be a string or number indicating the name or index of the desired statistic")
    }
  }

  index <- unique(drop(index))

  if (!several.ok && length(index) > 1L) {
    arg::err("{.arg {(.arg_index)}} must have length one")
  }

  if (is.numeric(index)) {
    if (!rlang::is_integerish(index)) {
      if (several.ok) {
        arg::err("{.arg {(.arg_index)}} must be a vector of positive integers")
      }
      else {
        arg::err("{.arg {(.arg_index)}} must be a positive integer")
      }
    }

    if (any(index > ncol(t))) {
      arg::err(sprintf("{.arg {(.arg_index)}} must be between 1 and %s", ncol(t)))
    }
  }
  else {
    if (is_null(colnames(t))) {
      arg::err("the estimates don't have names, so {.arg {(.arg_index)}} must be numeric")
    }

    index <- match(index, colnames(t))

    if (anyNA(index)) {
      arg::err(c("All entries in {.arg {(.arg_index)}} must be the names of available statistics to compute.",
                 " " = "The following are allowed: {.val {colnames(t)}}"))
    }
  }

  index
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

is_null <- function(x) {isTRUE(length(x) == 0L)}
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

#Draw `R` independent L'Ecuyer-CMRG streams, one per bootstrap replicate.
#
#Seeding each replicate's weights from its own stream is what makes the weights a
#function of the replicate index alone, rather than of how a backend happened to
#chunk and seed the work. That is what lets `fwb.array()` recover them without
#re-running anything, and what makes the results identical across every value of
#`cl` and `verbose`.
#
#Exactly one draw is taken from the caller's stream, so successive calls to `fwb()`
#without an intervening `set.seed()` use different streams, and the caller's RNG kind
#is never changed -- the L'Ecuyer chain is built inside `with_seed_preserved()`.
make_stream_seeds <- function(R) {

  if (identical(RNGkind()[1L], "L'Ecuyer-CMRG")) {
    #Consume a draw so repeated calls differ, then start the chain from the state
    #that draw left behind.
    sample.int(.Machine$integer.max, 1L)
    seed <- get(".Random.seed", globalenv(), mode = "integer", inherits = FALSE)
  }
  else {
    s <- sample.int(.Machine$integer.max, 1L)

    with_seed_preserved({
      suppressWarnings(RNGkind("L'Ecuyer-CMRG"))
      set.seed(s)
      seed <- get(".Random.seed", globalenv(), mode = "integer", inherits = FALSE)
    })
  }

  seeds <- matrix(0L, nrow = R, ncol = length(seed))

  for (i in seq_len(R)) {
    seed <- parallel::nextRNGStream(seed)
    seeds[i, ] <- seed
  }

  seeds
}

#Make the package's exported functions findable on `cluster` workers.
#
#`library()` is used rather than `requireNamespace()` because the lookup that fails
#happens through the search path: `statistic` is a closure whose enclosing environment
#is the worker's (empty) global environment, so loading the namespace without attaching
#it does not help. Errors are swallowed because a worker that cannot see the package
#library is only a problem for a `statistic` that needs it, and that failure reports
#itself clearly.
attach_on_workers <- function(cl) {
  try(parallel::clusterCall(cl, function() {
    suppressWarnings(suppressMessages(library("fwb", character.only = TRUE)))
    NULL
  }), silent = TRUE)

  invisible(NULL)
}

#Install a stream recorded by `make_stream_seeds()`. Assigning `.Random.seed` also
#sets the RNG kind from its first element, which is why every caller runs inside
#`with_seed_preserved()`.
set_stream_seed <- function(seed) {
  assign(".Random.seed", value = seed, envir = globalenv())
}

#Make sure `.Random.seed` exists before anything tries to read it.
#
#A session has no `.Random.seed` until the RNG is first used, so reading it in a fresh
#`Rscript` process fails with "object '.Random.seed' not found" -- reachable whenever
#nothing has drawn a random number yet, including a `statistic` that never does. Any RNG
#use initializes it from the current kind.
#
#Every read of `.Random.seed` in the package must be preceded by this, either directly or
#by being inside `with_seed_preserved()`. Keeping the check in one place is deliberate: it
#was previously inlined in two, and one copy was lost in a refactor.
ensure_random_seed <- function() {
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    runif(1L)
  }

  invisible(NULL)
}

with_seed_preserved <- function(expr, new_seed = NULL) {

  ensure_random_seed()

  old_seed <- list(random_seed = get(".Random.seed", globalenv(), mode = "integer",
                                     inherits = FALSE),
                   rng_kind = RNGkind())

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
