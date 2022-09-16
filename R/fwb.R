#' Fractional Weighted Bootstrap
#'
#' `fwb()` implements the fractional (random) weighted bootstrap, also known as the Bayesian bootstrap. Rather than resampling units to include in bootstrap samples, weights are drawn to be applied to a weighted estimator.
#'
#' @param data the dataset used to compute the statistic
#' @param statistic a function, which, when applied to `data`, returns a vector containing the statistic(s) of interest. The function should take at least two arguments; the first argument should correspond to the dataset and the second argument should correspond to a vector of weights. Any further arguments can be passed to `statistic` through the `...` argument.
#' @param R the number of bootstrap replicates. Default is 999 but more is always better. For the percentile bootstrap confidence interval to be exact, it can be beneficial to use one less than a multiple of 100.
#' @param cluster optional; a vector containing cluster membership. If supplied, will run the cluster bootstrap. See Details. Evaluated first in `data` and then in the global environment.
#' @param simple `logical`; if `TRUE`, weights will be computed on-the-fly in each bootstrap replication rather than all at once. This can save memory at the cost of some speed.
#' @param verbose `logical`; whether to display a progress bar.
#' @param cl a cluster object created by \pkgfun{parallel}{makeCluster}, or an integer to indicate the number of child-processes (integer values are ignored on Windows) for parallel evaluations. See \pkgfun{pbapply}{pblapply} for details. If `NULL`, no parallelization will take place.
#' @param ... other arguments passed to `statistic`.
#'
#' @return A `fwb` object, which also inherits from `boot`, with the following components:
#' \item{t0}{The observed value of `statistic` applied to `data` with uniform weights.}
#' \item{t}{A matrix with `R` rows, each of which is a bootstrap replicate of the result of calling `statistic`.}
#' \item{R}{The value of `R` as passed to `fwb()`.}
#' \item{data}{The `data` as passed to `fwb()`.}
#' \item{seed}{The value of `.Random.seed` just prior to generating the weights (after the first call to `statistic` with uniform weights).}
#' \item{statistic}{The function `statistic` as passed to `fwb()`.}
#' \item{call}{The original call to `fwb()`.}
#' \item{cluster}{The vector passed to `cluster`, if any.}
#'
#' @details `fwb()` implements the fractional weighted bootstrap and is meant to function as a drop-in for `boot::boot(., stype = "f")` (i.e., the usual bootstrap but with frequency weights representing the number of times each unit is drawn). In each bootstrap replication, the weights are sampled from independent exponential distributions with rate parameter 1 and then normalized to have a mean of 1, equivalent to drawing the weights from a Dirichlet distribution. The function supplied to `statistic` must incorporate the weights to compute a weighted statistic. For example, if the output is a regression coefficient, the weights supplied to the `w` argument of `statistic` should be supplied to the `weights` argument of `lm()`. These weights should be used any time frequency weights would be, since they are meant to function like frequency weights (which, in the case of the traditional bootstrap, would be integers). Unfortunately, there is no way for `fwb()` to know whether you are using the weights correctly, so care should be taken to ensure weights are correctly incorporated into the estimator.
#'
#' When fitting binomial regression models (e.g., logistic) using [glm()], it may be useful to change the `family` to a "quasi" variety (e.g., [quasibinomial()]) to avoid a spurious warning about "non-integer #successes".
#'
#' The cluster bootstrap can be requested by supplying a vector of cluster membership to `cluster`. Rather than generating a weight for each unit, a weight is generated for each cluster and then applied to all units in that cluster.
#'
#' Ideally, `statistic` should not involve a random element, or else it will not be straightforward to replicate the bootstrap results using the `seed` included in the output object. Setting a seed using [set.seed()] is always advised.
#'
#' The `print()` method displays the value of the statistics, the bias (the difference between the statistic and the mean of its bootstrap distribution), and the standard error (the standard deviation of the bootstrap distribution).
#'
#' @seealso [fwb.ci()] for calculating confidence intervals; [summary.fwb()] for displaying output in a clean way; [plot.fwb()] for plotting the bootstrap distributions; [vcovFWB()] for estimating the covariance matrix of estimates using the FWB; \pkgfun{boot}{boot} for the traditional bootstrap.
#'
#' @export
#'
#' @examplesIf requireNamespace("survival", quietly = TRUE)
#' # Performing a Weibull analysis of the Bearing Cage
#' # failure data as done in Xu et al. (2020)
#' data("bearingcage")
#'
#' weibull_est <- function(data, w) {
#'   fit <- survival::survreg(survival::Surv(hours, failure) ~ 1,
#'                            data = data, weights = w,
#'                            dist = "weibull")
#'
#'   c(eta = unname(exp(coef(fit))), beta = 1/fit$scale)
#' }
#'
#' boot_est <- fwb(bearingcage, statistic = weibull_est,
#'                 R = 199, verbose = FALSE)
#' boot_est
#'
#' #Get standard errors and CIs; uses bias-corrected
#' #percentile CI by default
#' summary(boot_est, ci.type = "bc")
#'
#' #Plot statistic distributions
#' plot(boot_est, index = "beta", type = "hist")
#'
#'
fwb <- function(data, statistic, R = 999, cluster = NULL, simple = FALSE, verbose = TRUE, cl = NULL, ...) {

  bcall <- match.call()

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
  # if (!all(c("data", "w") %in% names(formals(statistic)))) {
  #   chk::err("`statistic` must have a `data` argument and a `w` argument")
  # }

  #Check R
  chk::chk_count(R)

  clus <- substitute(cluster)
  cluster <- eval(clus, data, parent.frame())
  chk::chk_atomic(cluster)

  #Check simple
  chk::chk_flag(simple)

  chk::chk_flag(verbose)

  n <- nrow(data)
  if (is.null(n) || !chk::vld_count(n)) {
    chk::err("`data` must be present")
  }

  #Test fun
  t0 <- try(statistic(data, rep(1, n), ...))
  if (inherits(t0, "try-error")) {
    chk::err("there was an error running the function supplied to `statistic` on unit-weighted data. Error produced:\n\t",
             conditionMessage(attr(t0, "condition")))
  }
  if (!is.numeric(t0) || !is.null(dim(t0))) {
    chk::err("the output of the function supplied to `statistic` must be a numeric vector")
  }

  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE))
    runif(1)
  seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)

  if (is.null(cluster)) {
    if (simple) {
      FUN <- function(i) {
        w <- rexp(n)
        w <- n*w/sum(w)
        statistic(data, w, ...)
      }
    }
    else {
      w <- matrix(rexp(n * R), nrow = R, ncol = n, byrow = TRUE)
      w <- n*w/rowSums(w)
      FUN <- function(i) {
        statistic(data, w[i,], ...)
      }
    }
  }
  else {
    cluster <- factor(cluster)
    nc <- nlevels(cluster)
    cluster_numeric <- as.integer(cluster)

    if (simple) {
      FUN <- function(i) {
        cluster_w <- rexp(nc)
        cluster_w <- nc * cluster_w / sum(cluster_w)
        w <- cluster_w[cluster_numeric]
        statistic(data, w, ...)
      }
    }
    else {
      cluster_w <- matrix(rexp(R * nc), nrow = R, ncol = nc, byrow = TRUE)
      cluster_w <- nc * cluster_w / rowSums(cluster_w)
      w <- matrix(0, nrow = R, ncol = n)
      for (k in seq_len(R)) {
        w[k,] <- cluster_w[k, cluster_numeric]
      }
      FUN <- function(i) {
        statistic(data, w[i,], ...)
      }
    }
  }

  opb <- pbapply::pboptions(type = if (verbose) "timer" else "none")
  on.exit(pbapply::pboptions(opb))

  #Run bootstrap
  t <- do.call("rbind", pbapply::pblapply(seq_len(R), FUN, cl = cl))

  out <- list(t0 = t0,
              t = t,
              R = R,
              data = data,
              seed = seed,
              statistic = statistic,
              call = bcall,
              cluster = cluster)

  class(out) <- c("fwb", "boot")

  out
}

#' @describeIn fwb Print an `fwb` object
#'
#' @param x an `fwb` object; the output of a call to `fwb()`.
#' @param digits the number of significant digits to print
#' @param index the index or indices of the position of the quantity of interest in `x$t0` if more than one was specified in `fwb()`. Default is to print all quantities.
#'
#' @export
print.fwb <- function(x, digits = getOption("digits"), index = 1L:ncol(x$t), ...) {
  cl <- x$call
  t <- matrix(x$t[, index], nrow = nrow(x$t))
  allNA <- apply(t, 2L, function(t) all(is.na(t)))
  ind1 <- index[allNA]
  index <- index[!allNA]
  t <- matrix(t[, !allNA], nrow = nrow(t))
  rn <- paste("t", index, "*", sep = "")

  if (length(index) == 0L)
    op <- NULL
  else if (is.null(t0 <- x$t0)) {
    op <- cbind(apply(t, 2L, mean, na.rm = TRUE),
                apply(t, 2L, sd, na.rm = TRUE))
    dimnames(op) <- list(rn, c("mean", "std. error"))
  }
  else {
    t0 <- x$t0[index]
    op <- cbind(t0,
                apply(t, 2L, mean, na.rm = TRUE) - t0,
                apply(t, 2L, sd, na.rm = TRUE))
    dimnames(op) <- list(rn, c("original", "bias", "std. error"))

  }

  cat("FRACTIONAL WEIGHTED BOOTSTRAP\n")

  cat("\nCall:\n")
  dput(cl, control = NULL)
  cat("\nBootstrap Statistics :\n")
  if (!is.null(op))
    print(op, digits = digits)
  if (length(ind1) > 0L)
    for (j in ind1) cat(sprintf("WARNING: All values of t%s* are NA\n", j))
  invisible(x)
}

