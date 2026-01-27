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
    if (!((is.numeric(index) && chk::vld_whole_number(index) && index == 1) ||
          (chk::vld_string(index) && !is_null(colnames(t)) && index == colnames(t)[1L]))) {
      .wrn("only one statistic is available; ignoring {.arg index}")
    }

    return(1L)
  }

  if (!((is.character(index) || is.numeric(index)) && is_null(dim(index)))) {
    if (several.ok) {
      .err("{.arg index} must be a character or numeric vector indicating the names or indices of the desired statistics")
    }
    else {
      .err("{.arg index} must be a string or number indicating the name or index of the desired statistic")
    }
  }

  index <- unique(drop(index))

  if (!several.ok && length(index) > 1L) {
    .err("{.arg index} must have length one")
  }

  if (is.numeric(index)) {
    if (!chk::vld_whole_numeric(index)) {
      if (several.ok) {
        .err("{.arg index} must be a vector of positive integers")
      }
      else {
        .err("{.arg index} must be a positive integer")
      }
    }

    if (any(index > ncol(t))) {
      .err(sprintf("{.arg index} must be between 1 and %s", ncol(t)))
    }
  }
  else {
    if (is_null(colnames(t))) {
      .err("the estimates don't have names, so {.arg index} must be numeric")
    }

    index <- match(index, colnames(t))

    if (anyNA(index)) {
      .err("all entries in {.arg index} must be the names of available statistics to compute. The following are allowed: {.val {colnames(t)}}")
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

  if (chk::vld_string(quotes)) {
    return(paste0(quotes, x, str_rev(quotes)))
  }

  if (!chk::vld_count(quotes) || quotes > 2L) {
    stop("`quotes` must be boolean, 1, 2, or a string.")
  }

  if (quotes == 0L) {
    return(x)
  }

  x <- {
    if (quotes == 1L) sprintf("'%s'", x)
    else sprintf('"%s"', x)
  }

  x
}

str_rev <- function(x) {
  vapply(lapply(strsplit(x, NULL), rev), paste, character(1L), collapse = "")
}

#More informative and cleaner version of base::match.arg(). Uses chk and cli.
match_arg <- function(arg, choices, several.ok = FALSE, context = NULL) {
  #Replaces match.arg() but gives cleaner error message and processing
  #of arg.
  if (missing(arg)) {
    .err("no argument was supplied to {.fun match_arg} (this is a bug)")
  }

  arg.name <- deparse1(substitute(arg), width.cutoff = 500L)

  if (missing(choices)) {
    sysP <- sys.parent()
    formal.args <- formals(sys.function(sysP))
    choices <- eval(formal.args[[as.character(substitute(arg))]],
                    envir = sys.frame(sysP))
  }

  if (is_null(arg)) {
    return(choices[1L])
  }

  if (several.ok) {
    chk::chk_character(arg, x_name = add_quotes(arg.name, "`"))
  }
  else {
    chk::chk_string(arg, x_name = add_quotes(arg.name, "`"))

    if (identical(arg, choices)) {
      return(arg[1L])
    }
  }

  i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)

  if (all(i == 0L)) {
    one_of <- {
      if (length(choices) <= 1L) NULL
      else if (several.ok) "at least one of"
      else "one of"
    }

    .err("{context} the argument to {.arg {arg.name}} should be {one_of} {.or {.val {choices}}}")
  }

  i <- i[i > 0L]

  choices[i]
}

.tail <- function(x, n = 1L) {
  chk::chk_count(n)
  chk::chk_gt(n, 0)

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

pkg_caller_call <- function() {
  pn <- utils::packageName()
  package.funs <- c(getNamespaceExports(pn),
                    .getNamespaceInfo(asNamespace(pn), "S3methods")[, 3L])

  for (i in seq_len(sys.nframe())) {
    e <- sys.call(i)

    n <- rlang::call_name(e)

    if (is_not_null(n) && n %in% package.funs) {
      return(e)
    }
  }

  NULL
}

.err <- function(m, n = NULL, tidy = TRUE, cli = TRUE) {
  if (cli) {
    m <- eval.parent(substitute(cli::format_inline(.m), list(.m = m)))
  }

  chk::message_chk(m, n = n, tidy = tidy) |>
    cli::ansi_strwrap() |>
    paste(collapse = "\n") |>
    rlang::abort(call = pkg_caller_call())
}
.wrn <- function(m, n = NULL, tidy = TRUE, immediate = TRUE, cli = TRUE) {
  if (cli) {
    m <- eval.parent(substitute(cli::format_inline(.m), list(.m = m)))
  }

  m <- chk::message_chk(m, n = n, tidy = tidy)

  if (immediate && isTRUE(all.equal(0, getOption("warn")))) {
    rlang::with_options({
      m |>
        cli::ansi_strwrap() |>
        paste(collapse = "\n") |>
        rlang::warn()
    }, warn = 1)
  }
  else {
    m |>
      cli::ansi_strwrap() |>
      paste(collapse = "\n") |>
      rlang::warn()
  }
}
.msg <- function(m, n = NULL, tidy = TRUE, cli = TRUE) {
  if (cli) {
    m <- eval.parent(substitute(cli::format_inline(.m), list(.m = m)))
  }

  chk::message_chk(m, n = n, tidy = tidy) |>
    cli::ansi_strwrap() |>
    paste(collapse = "\n") |>
    rlang::inform(tidy = FALSE)
}

.chk_atomic_vector <- function(x, x_name = NULL) {
  if (is.atomic(x) && !is.matrix(x) && !is.array(x)) {
    return(invisible(x))
  }

  if (is_null(x_name)) {
    x_name <- chk::deparse_backtick_chk(substitute(x))
  }

  chk::abort_chk(x_name, " must be an atomic vector", x = x)
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
