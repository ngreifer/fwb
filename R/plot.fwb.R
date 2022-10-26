#' Plots of the Output of a Fractional Weighted Bootstrap
#'
#' `plot.fwb()` takes an `fwb` object and produces plots for the bootstrap replicates of the statistic of interest.
#'
#' @param x an `fwb` object; the output of a call to [fwb()].
#' @param index the index of the position of the quantity of interest in `x$t0` if more than one was specified in `fwb()`. Only one value is allowed at a time. By default the first statistic is used.
#' @param qdist `character`; when a Q-Q plot is requested (as it is by default; see `type` argument below), the distribution against which the Q-Q plot should be drawn. Allowable options include `"norm"` (normal distribution - the default) and `"chisq"` (chi-squared distribution).
#' @param nclass when a histogram is requested (as it is by default; see `type` argument below), the number of classes to be used. The default is the integer between 10 and 100 closest to `ceiling(length(R)/25)` where `R` is the number of bootstrap replicates.
#' @param df if `qdist` is `"chisq"`, the degrees of freedom for the chi-squared distribution to be used. If not supplied, the degrees of freedom will be estimated using maximum likelihood.
#' @param type the type of plot to display. Allowable options include `"hist"` for a histogram of the bootstrap estimates and `"qq"` for a Q-Q plot of the estimates against the distribution supplied to `qdist`.
#' @param ... ignored.
#'
#' @return `x` is returned invisibly.
#'
#' @details This function can produces two side-by-side plots: a histogram of the bootstrap replicates and a Q-Q plot of the bootstrap replicates against theoretical quantiles of a supplied distribution (normal or chi-squared). For the histogram, a vertical dotted line indicates the position of the estimate computed in the original sample. For the Q-Q plot, the expected line is plotted.
#'
#' @seealso [fwb()], [summary.fwb()], \pkgfun{boot}{plot.boot}, [hist()], [qqplot()]
#'
#' @export
#'
#' @examples
#' # See examples at help("fwb")
plot.fwb <- function(x, index = 1, qdist = "norm", nclass = NULL, df, type = c("hist", "qq"), ...) {

  index <- check_index(index, x[["t"]])

  chk::chk_string(qdist)

  t <- x[["t"]][, index]
  t0 <- x[["t0"]][index]

  t <- t[is.finite(t)]
  if (all_the_same(t)) {
    chk::wrn("All values of t* are equal to ", mean(t, na.rm = TRUE))
    return(invisible(x))
  }

  chk::chk_character(type)
  type <- tryCatch(match.arg(tolower(type), c("hist", "qq"), several.ok = TRUE),
                   error = function(e) {
                     .err("`type` must be one or more of \"hist\" or \"qq\"")
                   })

  opar <- graphics::par(mfrow = c(1, length(type)))
  on.exit(graphics::par(opar))

  for (i in type) {
    if (i == "hist") {
      if (is.null(nclass)) {
        nclass <- min(max(ceiling(length(t)/25), 10), 100)
      }
      else {
        chk::chk_count(nclass)
      }

      rg <- range(t)
      if (t0 < rg[1L])
        rg[1L] <- t0
      else if (t0 > rg[2L])
        rg[2L] <- t0
      rg <- rg + 0.05 * c(-1, 1) * diff(rg)
      lc <- diff(rg)/(nclass - 2)
      n1 <- ceiling((t0 - rg[1L])/lc)
      n2 <- ceiling((rg[2L] - t0)/lc)
      bks <- t0 + (-n1:n2) * lc

      graphics::hist(t, breaks = bks, probability = TRUE,
                     xlab = colnames(x[["t"]])[index],
                     main = sprintf("Histogram of %s", colnames(x[["t"]])[index]))

      graphics::abline(v = t0, lty = 2)
    }
    else if (i == "qq") {
      p <- ppoints(x$R, a = .5)

      if (qdist == "chisq") {
        if (missing(df)) {
          df <- estimate_chisq_df(t)
        }
        else {
          chk::chk_number(df)
        }
        qf <- function(p_) qchisq(p_, df = df)
        qlab <- sprintf("Quantiles of Chi-squared(%s)", round(df, 2))
      }
      else if (qdist == "norm") {
        qf <- function(p_) qnorm(p_)
        qlab <- "Quantiles of Standard Normal"
      }
      else {
        .err(sprintf("\"%s\" distribution not supported: using normal instead",
                         qdist))
      }

      qqplot(qf(p), t, xlab = qlab, ylab = colnames(x[["t"]])[index])
      qqline(t, distribution = qf, qtype = 5, lty = 2)
    }
  }

  invisible(x)
}

#Estimates df of chisq distribution using MLE; port of MASS::fitdsitr
estimate_chisq_df <- function(x) {
  myfn <- function(df, x, ...) -sum(log(dchisq(x, df, ...)))

  res <- optimize(myfn, interval = c(0, 10 * mean(x)), x = x,
                  tol = 1e-8)
  res$minimum
}
