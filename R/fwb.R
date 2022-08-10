fwb <- function(data, statistic, R = 999, simple = FALSE, cl = NULL, ...) {

  #Check data
  if (missing(data)) {
    chk::err("`data` must be specified")
  }
  chk::chk_data(data)

  #Check statistic
  if (missing(statistic)) {
    chk::err("`statistic` must be specified")
  }
  chk::chk_function(statistic)
  if (!all(c("data", "w") %in% names(formals(statistic)))) {
    chk::err("`statistic` must have a `data` argument and a `w` argument")
  }

  #Check R
  chk::chk_count(R)

  #Check simple
  chk::chk_flag(simple)

  n <- nrow(data)

  #Test fun
  t0 <- try(statistic(data = data, w = rep(1, n), ...))
  if (inherits(t0, "try-error")) {
    chk::err("there was an error running the function supplied to `statistic` on unit-weighted data. Error produced:\n\t",
             conditionMessage(attr(t0, "condition")))
  }
  if (!is.numeric(t0) || !is.null(dim(t0))) {
    chk::err("the output of the function supplied to `statistic` must be a numeric vector")
  }
##
  if (simple) {
    FUN <- function(i) {
        w <- rexp(n)
        w <- n*w/sum(w)
        statistic(data = data, w = w, ...)
      }
  }
  else {
    w <- matrix(rexp(n * R), nrow = n, ncol = R)
    w <- sweep(w, 2, colMeans(w), "/")
    FUN <- function(i) {
      statistic(data = data, w = w[,i], ...)
    }
  }


  ests <- do.call("rbind", pbapply::pblapply(seq_len(R), FUN, cl = cl))

  out <- list(ests = ests,
              t0 = statistic(rep(1, n)))

  class(out) <- "fwb"

  out
}

fwbsum <- function(f) {
  out <- matrix(nrow = ncol(f$ests), ncol = 4,
                dimnames = list(colnames(f$ests),
                                c("Estimate", "Std. Error", "CI 2.5 %", "CI 97.5 %")))

  out[,"Estimate"] <- f$t0
  out[,"Std. Error"] <- apply(f$ests, 2, sd)
  out[,"CI 2.5 %"] <- apply(f$ests, 2, quantile, probs = .025)
  out[,"CI 97.5 %"] <- apply(f$ests, 2, quantile, probs = .975)

  class(out) <- c("fwbsum", class(out))
  out
}


