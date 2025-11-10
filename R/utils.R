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
      chk::wrn("only one statistic is available; ignoring `index`")
    }

    return(1L)
  }

  if (!((is.character(index) || is.numeric(index)) && is_null(dim(index)))) {
    if (several.ok) {
      .err("`index` must be a character or numeric vector indicating the names or indices of the desired statistics")
    }
    else {
      .err("`index` must be a string or number indicating the name or index of the desired statistic")
    }
  }

  index <- unique(drop(index))

  if (!several.ok && length(index) > 1L) {
    .err("`index` must have length one")
  }

  if (is.numeric(index)) {
    if (!chk::vld_whole_numeric(index)) {
      if (several.ok) {
        .err("`index` must be a vector of positive integers")
      }
      else {
        .err("`index` must be a positive integer")
      }
    }

    if (any(index > ncol(t))) {
      .err(sprintf("`index` must be between 1 and %s", ncol(t)))
    }
  }
  else {
    if (is_null(colnames(t))) {
      .err("the estimates don't have names, so `index` must be numeric")
    }

    index <- match(index, colnames(t))

    if (anyNA(index)) {
      .err("all entries in `index` must be the names of available statistics to compute.\n  The following are allowed: ",
           toString(add_quotes(colnames(t), 2L)))
    }
  }

  index
}

word_list <- function(word.list = NULL, and.or = "and", is.are = FALSE, quotes = FALSE) {
  #When given a vector of strings, creates a string of the form "a and b"
  #or "a, b, and c"
  #If is.are, adds "is" or "are" appropriately

  word.list <- setdiff(word.list, c(NA_character_, ""))

  if (is_null(word.list)) {
    out <- ""
    attr(out, "plural") <- FALSE
    return(out)
  }

  word.list <- add_quotes(word.list, quotes)

  L <- length(word.list)

  if (L == 1L) {
    out <- word.list
    if (is.are) out <- paste(out, "is")
    attr(out, "plural") <- FALSE
    return(out)
  }

  if (is_null(and.or) || isFALSE(and.or)) {
    out <- toString(word.list)
  }
  else {
    and.or <- match_arg(and.or, c("and", "or"))

    if (L == 2L) {
      out <- sprintf("%s %s %s",
                     word.list[1L],
                     and.or,
                     word.list[2L])
    }
    else {
      out <- sprintf("%s, %s %s",
                     toString(word.list[-L]),
                     and.or,
                     word.list[L])
    }
  }

  if (is.are) {
    out <- sprintf("%s are", out)
  }

  attr(out, "plural") <- TRUE

  out
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

match_arg <- function(arg, choices, several.ok = FALSE, context = NULL) {
  #Replaces match.arg() but gives cleaner error message and processing
  #of arg.
  if (missing(arg)) {
    stop("No argument was supplied to match_arg.")
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
    .err(sprintf("%sthe argument to `%s` should be %s%s",
                 if (is_null(context)) "" else sprintf("%s ", context),
                 arg.name,
                 ngettext(length(choices), "", if (several.ok) "at least one of " else "one of "),
                 word_list(choices, and.or = "or", quotes = 2L)))
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

  class(x) <- c(.class, class(x))

  x
}

pkg_caller_call <- function() {
  pn <- utils::packageName()
  package.funs <- c(getNamespaceExports(pn),
                    .getNamespaceInfo(asNamespace(pn), "S3methods")[, 3L])

  for (i in seq_len(sys.nframe())) {
    e <- sys.call(i)

    n <- rlang::call_name(e)

    if (is_null(n)) {
      next
    }

    if (n %in% package.funs) {
      return(e)
    }
  }

  NULL
}

.err <- function(..., n = NULL, tidy = TRUE) {
  m <- chk::message_chk(..., n = n, tidy = tidy)
  rlang::abort(paste(strwrap(m), collapse = "\n"),
               call = pkg_caller_call())
}
.wrn <- function(..., n = NULL, tidy = TRUE, immediate = TRUE) {
  m <- chk::message_chk(..., n = n, tidy = tidy)

  if (immediate && isTRUE(all.equal(0, getOption("warn")))) {
    rlang::with_options({
      rlang::warn(paste(strwrap(m), collapse = "\n"))
    }, warn = 1)
  }
  else {
    rlang::warn(paste(strwrap(m), collapse = "\n"))
  }
}
.msg <- function(..., n = NULL, tidy = TRUE) {
  m <- chk::message_chk(..., n = n, tidy = tidy)
  rlang::inform(paste(strwrap(m), collapse = "\n"), tidy = FALSE)
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
