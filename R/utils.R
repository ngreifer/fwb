make_gen_weights <- function(wtype) {
  wtype <- tolower(wtype)
  wtype <- match_arg(wtype, c("exp", "multinom", "poisson", "mammen"))

  fun <- switch(wtype,
         "exp" = function(n, R) {
           w <- matrix(rexp(n * R), nrow = R, ncol = n, byrow = TRUE)
           n * w/rowSums(w)
         },
         "poisson" = function(n, R) {
           matrix(rpois(n * R, 1), nrow = R, ncol = n, byrow = TRUE)
         },
         # "multinom" = function(n, R) {
         #   rmultinom(R, n, rep(1/n, n))
         # },
         "multinom" = function(n, R) {
           i <- sample.int(n, n * R, replace = TRUE)
           dim(i) <- c(R, n)
           t(apply(i, 1, tabulate, n))
         },
         "mammen" = function(n, R) {
           sqrt5 <- sqrt(5)
           w <- matrix((3-sqrt5)/2 + sqrt5 * rbinom(n * R, 1, .5 - 1/(2*sqrt5)),
                       nrow = R, ncol = n, byrow = TRUE)
           n * w/rowSums(w)
         })

  attr(fun, "wtype") <- wtype
  fun
}

na.rem <- function(x) {
  #A faster na.omit for vectors
  x[!is.na(x)]
}

check_if_zero <- function(x) {
  # this is the default tolerance used in all.equal
  tolerance <- .Machine$double.eps^0.5
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
  else all(x == x[1])
}

#Format percentage for CI labels
fmt.prc <- function(probs, digits = 3) {
  paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), "%")
}

check_index <- function(index, t, several.ok = FALSE) {
  if (is_null(index)) {
    return(1L)
  }

  if (ncol(t) == 1L) {
    if (!((is.numeric(index) && chk::vld_whole_number(index) && index == 1) ||
          (chk::vld_string(index) && !is.null(colnames(t)) && index == colnames(t)[1]))) {
      chk::wrn("only one statistic is available; ignoring `index`")
    }

    return(1L)
  }

  if (!((is.character(index) || is.numeric(index)) && is_null(dim(index)))) {
    if (!several.ok) {
      .err("`index` must be a string or number indicating the name or index of the desired statistic")
    }
    else {
      .err("`index` must be a character or numeric vector indicating the names or indices of the desired statistics")
    }
  }

  index <- unique(drop(index))

  if (!several.ok && length(index) > 1) {
    .err("`index` must have length one")
  }

  if (is.numeric(index)) {
    if (!chk::vld_whole_numeric(index)) {
      if (!several.ok) {
        .err("`index` must be a positive integer")
      }
      else {
        .err("`index` must be a vector of positive integers")
      }
    }

    if (any(index > ncol(t))) {
      .err("`index` must be between 1 and ", ncol(t))
    }
  }
  else {
    if (is.null(colnames(t))) {
      .err("the estimates don't have names, so `index` must be numeric")
    }

    index <- match(index, colnames(t))

    if (anyNA(index)) {
      .err("all entries in `index` must be the names of available statistics to compute.\n  The following are allowed: ", paste(add_quotes(colnames(t), 2L), collapse = ", "))
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
    out <- paste(word.list, collapse = ", ")
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
                     paste(word.list[-L], collapse = ", "),
                     and.or,
                     word.list[L])
    }
  }

  if (is.are) out <- sprintf("%s are", out)

  attr(out, "plural") <- TRUE

  out
}

add_quotes <- function(x, quotes = 2L) {
  if (isFALSE(quotes)) {
    return(x)
  }

  if (isTRUE(quotes))
    quotes <- '"'

  if (chk::vld_string(quotes)) {
    return(paste0(quotes, x, quotes))
  }

  if (!chk::vld_count(quotes) || quotes > 2) {
    stop("`quotes` must be boolean, 1, 2, or a string.")
  }

  if (quotes == 0L) {
    return(x)
  }

  x <- {
    if (quotes == 1) sprintf("'%s'", x)
    else sprintf('"%s"', x)
  }

  x
}

match_arg <- function(arg, choices, several.ok = FALSE) {
  #Replaces match.arg() but gives cleaner error message and processing
  #of arg.
  if (missing(arg)) {
    stop("No argument was supplied to match_arg.")
  }

  arg.name <- deparse1(substitute(arg), width.cutoff = 500L)

  if (missing(choices)) {
    formal.args <- formals(sys.function(sysP <- sys.parent()))
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
  if (all(i == 0L))
    .err(sprintf("the argument to `%s` should be %s%s",
                 arg.name,
                 ngettext(length(choices), "", if (several.ok) "at least one of " else "one of "),
                 word_list(choices, and.or = "or", quotes = 2)))
  i <- i[i > 0L]

  choices[i]
}

is_null <- function(x) {length(x) == 0L}
is_not_null <- function(x) {!is_null(x)}

pkg_caller_call <- function() {
  pn <- utils::packageName()
  package.funs <- c(getNamespaceExports(pn),
                    .getNamespaceInfo(asNamespace(pn), "S3methods")[, 3])

  for (i in seq_len(sys.nframe())) {
    e <- sys.call(i)

    if (is_null(n <- rlang::call_name(e))) {
      next
    }

    if (n %in% package.funs) {
      return(e)
    }
  }

  NULL
}

#chk utilities
.err <- function(..., n = NULL, tidy = TRUE) {
  m <- chk::message_chk(..., n = n, tidy = tidy)
  rlang::abort(paste(strwrap(m), collapse = "\n"),
               call = pkg_caller_call())
}
.wrn <- function(..., n = NULL, tidy = TRUE, immediate = TRUE) {
  if (immediate && isTRUE(all.equal(0, getOption("warn")))) {
    op <- options(warn = 1)
    on.exit(options(op))
  }
  m <- chk::message_chk(..., n = n, tidy = tidy)
  rlang::warn(paste(strwrap(m), collapse = "\n"))
}
.msg <- function(..., n = NULL, tidy = TRUE) {
  m <- chk::message_chk(..., n = n, tidy = tidy)
  rlang::inform(paste(strwrap(m), collapse = "\n"), tidy = FALSE)
}

.chk_atomic_vector <- function(x, x_name = NULL) {
  if (is.atomic(x) && !is.matrix(x) && !is.array(x)) {
    return(invisible(x))
  }
  if (is.null(x_name))
    x_name <- chk::deparse_backtick_chk(substitute(x))
  chk::abort_chk(x_name, " must be an atomic vector", x = x)
}
