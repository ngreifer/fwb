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

#Check if all values are the same
all_the_same <- function(x) {
  if (is.numeric(x)) return(abs(max(x) - min(x)) < 1e-9)

  length(unique(x)) == 1
}

#Format percentage for CI labels
fmt.prc <- function(probs, digits = 3) {
  paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), "%")
}

check_index <- function(index, t, several.ok = FALSE) {
  if (length(index) == 0) return(1L)

  if (ncol(t) == 1) {
    if (!((is.numeric(index) && chk::vld_whole_number(index) && index == 1) ||
          (chk::vld_string(index) && !is.null(colnames(t)) && index == colnames(t)[1]))) {
      chk::wrn("only one statistic is available; ignoring `index`")
    }
    return(1L)
  }

  if (!((is.character(index) || is.numeric(index)) && is.null(dim(index)))) {
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
    if (!all(index %in% colnames(t))) {
      .err("all entries in `index` must be the names of available statistics to compute.\n  The following are allowed: ", paste(dQuote(colnames(t), FALSE), collapse = ", "))
    }
    index <- match(index, colnames(t))
  }

  index
}

add_quotes <- function(x, quotes = 2) {
  if (!isFALSE(quotes)) {
    if (isTRUE(quotes) || as.integer(quotes) == 2) x <- paste0("\"", x, "\"")
    else if (as.integer(quotes) == 1) x <- paste0("\'", x, "\'")
    else stop("'quotes' must be boolean, 1, or 2.")
  }
  x
}

word_list <- function(word.list = NULL, and.or = c("and", "or"), is.are = FALSE, quotes = FALSE) {
  #When given a vector of strings, creates a string of the form "a and b"
  #or "a, b, and c"
  #If is.are, adds "is" or "are" appropriately
  L <- length(word.list)
  word.list <- add_quotes(word.list, quotes)

  if (L == 0) {
    out <- ""
    attr(out, "plural") <- FALSE
  }
  else {
    word.list <- word.list[!word.list %in% c(NA_character_, "")]
    L <- length(word.list)
    if (L == 0) {
      out <- ""
      attr(out, "plural") <- FALSE
    }
    else if (L == 1) {
      out <- word.list
      if (is.are) out <- paste(out, "is")
      attr(out, "plural") <- FALSE
    }
    else {
      and.or <- match_arg(and.or)
      if (L == 2) {
        out <- paste(word.list, collapse = paste0(" ", and.or," "))
      }
      else {
        out <- paste(paste(word.list[seq_len(L-1)], collapse = ", "),
                     word.list[L], sep = paste0(", ", and.or," "))

      }
      if (is.are) out <- paste(out, "are")
      attr(out, "plural") <- TRUE
    }

  }
  return(out)
}

match_arg <- function(arg, choices, several.ok = FALSE) {
  #Replaces match.arg() but gives cleaner error message and processing
  #of arg.
  if (missing(arg))
    .err("No argument was supplied to `match_arg()`.", call. = FALSE)
  arg.name <- deparse1(substitute(arg), width.cutoff = 500L)

  if (missing(choices)) {
    formal.args <- formals(sys.function(sysP <- sys.parent()))
    choices <- eval(formal.args[[as.character(substitute(arg))]],
                    envir = sys.frame(sysP))
  }

  if (is.null(arg))
    return(choices[1L])

  if (!is.character(arg))
    .err(sprintf("The argument to `%s` must be `NULL` or a character vector", arg.name), call. = FALSE)

  if (!several.ok) {
    if (identical(arg, choices))
      return(arg[1L])
    if (length(arg) > 1L)
      .err(sprintf("The argument to `%s` must be of length 1", arg.name), call. = FALSE)
  }
  else if (length(arg) == 0)
    .err(sprintf("The argument to `%s` must be of length >= 1", arg.name), call. = FALSE)

  i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
  if (all(i == 0L))
    .err(sprintf("The argument to `%s` should be %s %s.",
                 arg.name, ngettext(length(choices), "", if (several.ok) "at least one of " else "one of "),
                 word_list(choices, and.or = "or", quotes = 2)),
         call. = FALSE)
  i <- i[i > 0L]
  if (!several.ok && length(i) > 1)
    .err("There is more than one match in `match_arg`")
  choices[i]
}

pkg_caller_call <- function(start = 1) {
  package.funs <- c(getNamespaceExports(utils::packageName()),
                    .getNamespaceInfo(asNamespace(utils::packageName()), "S3methods")[,3])
  k <- start #skip checking pkg_caller_call()
  e_max <- NULL
  while(!is.null(e <- rlang::caller_call(k))) {
    if (!is.null(n <- rlang::call_name(e)) &&
        n %in% package.funs) e_max <- k
    k <- k + 1
  }
  rlang::caller_call(e_max)
}

.err <- function(...) {
  chk::err(..., call = pkg_caller_call(start = 2))
}

.wrn <- function(..., immediate = TRUE) {
  if (immediate && isTRUE(all.equal(getOption("warn"), 0))) {
    op <- options(warn = 1)
    on.exit(options(op))
  }
  chk::wrn(...)
}

.chk_atomic_vector <- function(x, x_name = NULL) {
  if (is.atomic(x) && !is.matrix(x) && !is.array(x)) {
    return(invisible(x))
  }
  if (is.null(x_name))
    x_name <- chk::deparse_backtick_chk(substitute(x))
  chk::abort_chk(x_name, " must be an atomic vector", x = x)
}
