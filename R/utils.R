#Check if all values are the same
all_the_same <- function(x) {
  if (is.numeric(x)) return(abs(max(x) - min(x)) < 1e-9)
  return(length(unique(x)) == 1)
}

#Format percentage for CI labels
fmt.prc <- function(probs, digits = 3) {
  paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), "%")
}

check_index <- function(index, names, several.ok = FALSE) {
  if (length(index) == 0) return(1L)
  else if (length(names) == 1) {
    if (!((is.numeric(index) && chk::vld_whole_number(index) && index == 1) ||
          (chk::vld_string(index) && index == names[1]))) {
      chk::wrn("only one statistic is available; ignoring `index`")
    }
    return(1L)
  }

  if (!((is.character(index) || is.numeric(index)) && is.null(dim(index)))) {
    if (!several.ok) {
      chk::err("`index` must be a string or number indicating the name or index of the desired statistic")
    }
    else {
      chk::err("`index` must be a character or numeric vector indicating the names or indices of the desired statistics")
    }
  }

  index <- unique(drop(index))

  if (!several.ok && length(index) > 1) {
    chk::err("`index` must have length one")
  }

  if (is.numeric(index)) {
    if (!chk::vld_whole_numeric(index)) {
      if (!several.ok) {
        chk::err("`index` must be a positive integer")
      }
      else {
        chk::err("`index` must be a vector of positive integers")
      }
    }
    if (any(index > length(names))) {
      chk::err("`index` must be between 1 and ", length(names))
    }
  }
  else {
    if (!all(index %in% names)) {
      chk::err("all entries in `index` must be the names of available statistics to compute. The following are allowed:\n", paste(dQuote(names, FALSE)))
    }
    index <- match(index, names)
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
    chk::err("No argument was supplied to `match_arg()`.", call. = FALSE)
  arg.name <- deparse1(substitute(arg), width.cutoff = 500L)

  if (missing(choices)) {
    formal.args <- formals(sys.function(sysP <- sys.parent()))
    choices <- eval(formal.args[[as.character(substitute(arg))]],
                    envir = sys.frame(sysP))
  }

  if (is.null(arg))
    return(choices[1L])
  else if (!is.character(arg))
    chk::err(sprintf("The argument to `%s` must be `NULL` or a character vector", arg.name), call. = FALSE)
  if (!several.ok) {
    if (identical(arg, choices))
      return(arg[1L])
    if (length(arg) > 1L)
      chk::err(sprintf("The argument to `%s` must be of length 1", arg.name), call. = FALSE)
  }
  else if (length(arg) == 0)
    chk::err(sprintf("The argument to `%s` must be of length >= 1", arg.name), call. = FALSE)

  i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
  if (all(i == 0L))
    chk::err(sprintf("The argument to `%s` should be %s %s.",
                 arg.name, ngettext(length(choices), "", if (several.ok) "at least one of " else "one of "),
                 word_list(choices, and.or = "or", quotes = 2)),
         call. = FALSE)
  i <- i[i > 0L]
  if (!several.ok && length(i) > 1)
    chk::err("There is more than one match in `match_arg`")
  choices[i]
}
