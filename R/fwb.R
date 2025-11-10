#' Fractional Weighted Bootstrap
#'
#' `fwb()` implements the fractional (random) weighted bootstrap, also known as the Bayesian bootstrap. Rather than resampling units to include in bootstrap samples, weights are drawn to be applied to a weighted estimator.
#'
#' @param data the dataset used to compute the statistic
#' @param statistic a function, which, when applied to `data`, returns a vector containing the statistic(s) of interest. The function should take at least two arguments; the first argument should correspond to the dataset and the second argument should correspond to a vector of weights. Any further arguments can be passed to `statistic` through the `...` argument.
#' @param R the number of bootstrap replicates. Default is 999 but more is always better. For the percentile bootstrap confidence interval to be exact, it can be beneficial to use one less than a multiple of 100.
#' @param cluster optional; a vector containing cluster membership. If supplied, will run the cluster bootstrap. See Details. Evaluated first in `data` and then in the global environment.
#' @param simple `logical`; if `TRUE`, weights will be generated on-the-fly in each bootstrap replication; if `FALSE`, all weights will be generated at once and then supplied to `statistic`. The default (`NULL`) sets to `FALSE` if `wtype = "multinom"` and to `TRUE` otherwise.
#' @param wtype string; the type of weights to use. Allowable options include `"exp"` (the default), `"pois"`, `"multinom"`, `"mammen"`, `"beta"`, and `"power"`. See Details. See [set_fwb_wtype()] to set a global default.
#' @param strata optional; a vector containing stratum membership for stratified bootstrapping. If supplied, will essentially perform a separate bootstrap within each level of `strata`. This does not affect results when `wtype = "poisson"`.
#' @param drop0 `logical`; when `wtype` is `"multinom"` or `"poisson"`, whether to drop units that are given weights of 0 from the dataset and weights supplied to `statistic` in each iteration. If `NA`, weights of 0 will be set to `NA` instead. Ignored for other `wtype`s because they don't produce 0 weights. Default is `FALSE`.
#' @param verbose `logical`; whether to display a progress bar.
#' @param cl a cluster object created by \pkgfun{parallel}{makeCluster}, an integer to indicate the number of child-processes (integer values are ignored on Windows) for parallel evaluations, or the string `"future"` to use a `future` backend. See the `cl` argument of \pkgfun{pbapply}{pblapply} for details. If `NULL`, no parallelization will take place. See `vignette("fwb-rep")` for details.
#' @param ... other arguments passed to `statistic`.
#'
#' @returns
#' An `fwb` object, which also inherits from `boot`, with the following components:
#'
#' \item{t0}{The observed value of `statistic` applied to `data` with uniform weights.}
#' \item{t}{A matrix with `R` rows, each of which is a bootstrap replicate of the result of calling `statistic`.}
#' \item{R}{The value of `R` as passed to `fwb()`.}
#' \item{data}{The `data` as passed to `fwb()`.}
#' \item{seed}{The value of `.Random.seed` just prior to generating the weights (after the first call to `statistic` with uniform weights).}
#' \item{statistic}{The function `statistic` as passed to `fwb()`.}
#' \item{call}{The original call to `fwb()`.}
#' \item{cluster}{The vector passed to `cluster`, if any.}
#' \item{strata}{The vector passed to `strata`, if any.}
#' \item{wtype}{The type of weights used as determined by the `wtype` argument.}
#'
#' `fwb` objects have [coef()] and [vcov()] methods, which extract the `t0` component and covariance of the `t` components, respectively.
#'
#' @details `fwb()` implements the fractional weighted bootstrap and is meant to function as a drop-in for `boot::boot(., stype = "f")` (i.e., the usual bootstrap but with frequency weights representing the number of times each unit is drawn). In each bootstrap replication, when `wtype = "exp"` (the default), the weights are sampled from independent exponential distributions with rate parameter 1 and then normalized to have a mean of 1, equivalent to drawing the weights from a Dirichlet distribution. Other weights are allowed as determined by the `wtype` argument (see below for details). The function supplied to `statistic` must incorporate the weights to compute a weighted statistic. For example, if the output is a regression coefficient, the weights supplied to the `w` argument of `statistic` should be supplied to the `weights` argument of `lm()`. These weights should be used any time frequency weights would be, since they are meant to function like frequency weights (which, in the case of the traditional bootstrap, would be integers). Unfortunately, there is no way for `fwb()` to know whether you are using the weights correctly, so care should be taken to ensure weights are correctly incorporated into the estimator.
#'
#' When fitting binomial regression models (e.g., logistic) using [glm()], it may be useful to change the `family` to a "quasi" variety (e.g., [quasibinomial()]) to avoid a spurious warning about "non-integer #successes".
#'
#' The cluster bootstrap can be requested by supplying a vector of cluster membership to `cluster`. Rather than generating a weight for each unit, a weight is generated for each cluster and then applied to all units in that cluster.
#'
#' Bootstrapping can be performed within strata by supplying a vector of stratum membership to `strata`. This essentially rescales the weights within each stratum to have a mean of 1, ensuring that the sum of weights in each stratum is equal to the stratum size. For multinomial weights, using strata is equivalent to drawing samples with replacement from each stratum. Strata do not affect bootstrapping when using Poisson weights.
#'
#' Ideally, `statistic` should not involve a random element, or else it will not be straightforward to replicate the bootstrap results using the `seed` included in the output object. Setting a seed using [set.seed()] is always advised. See `vignette("fwb-rep")` for details.
#'
#' The `print()` method displays the value of the statistics, the bias (the difference between the statistic and the mean of its bootstrap distribution), and the standard error (the standard deviation of the bootstrap distribution).
#'
#' ## Weight types
#'
#' Different types of weights can be supplied to the `wtype` argument. A global default can be set using [set_fwb_wtype()]. The allowable weight types are described below.
#'
#' \describe{
#' \item{`"exp"`}{
#' Draws weights from an exponential distribution with rate parameter 1 using [rexp()]. These weights are the usual "Bayesian bootstrap" weights described in Xu et al. (2020). They are equivalent to drawing weights from a uniform Dirichlet distribution, which is what gives these weights the interpretation of a Bayesian prior. These positive weights have a mean and variance of 1 and skewness of 2. The weights are scaled to have a mean of 1 within each stratum (or in the full sample if `strata` is not supplied).
#' }
#' \item{`"multinom"`}{
#' Draws integer weights using [sample()], which samples unit indices with replacement and uses the tabulation of the indices as frequency weights. This is equivalent to drawing weights from a multinomial distribution. Using `wtype = "multinom"` is the same as using `boot::boot(., stype = "f")` instead of `fwb()` (i.e., the resulting estimates will be identical). When `strata` is supplied, unit indices are drawn with replacement within each stratum so that the sum of the weights in each stratum is equal to the stratum size.
#' }
#' \item{`"poisson"`}{
#' Draws integer weights from a Poisson distribution with 1 degree of freedom using [rpois()]. This is an alternative to the multinomial weights that yields similar estimates (especially as the sample size grows) but can be faster. Note `strata` is ignored when using `"poisson"`.
#' }
#' \item{`"mammen"`}{
#' Draws weights from a modification of the distribution described by Mammen (1983) for use in the wild bootstrap. These positive weights have a mean, variance, and skewness of 1, making them second-order accurate (in contrast to the usual exponential weights, which are only first-order accurate). The weights \eqn{w} are drawn such that \eqn{P(w=(3+\sqrt{5})/2)=(\sqrt{5}-1)/2\sqrt{5}} and \eqn{P(w=(3-\sqrt{5})/2)=(\sqrt{5}+1)/2\sqrt{5}}. The weights are scaled to have a mean of 1 within each stratum (or in the full sample if `strata` is not supplied).
#' }
#' \item{`"beta"`}{
#' Draws weights from a \eqn{\text{Beta}(1/2, 3/2)} distribution using [rbeta()] as described by Owen (2025). These positive weights have a mean, variance, and skewness of 1 when scaled by a factor of 4, making them second-order accurate. The weights are scaled to have a mean of 1 within each stratum (or in the full sample if `strata` is not supplied).
#' }
#' \item{`"power"`}{
#' Draws weights from a \eqn{\text{Beta}(\sqrt{2} - 1, 1)} distribution using [rbeta()] as described by Owen (2025). These positive weights have a mean and variance of 1 and skewness of \eqn{2(\sqrt{2} - 1)} when scaled by a factor of \eqn{2+\sqrt{2}}. The weights are scaled to have a mean of 1 within each stratum (or in the full sample if `strata` is not supplied).
#' }
#' }
#'
#' `"exp"` is the default due to it being the formulation described in Xu et al. (2020) and in the most formulations of the Bayesian bootstrap; it should be used if one wants to remain in line with these guidelines or to maintain a Bayesian flavor to the analysis, whereas `"mammen"`, `"beta"`, and `"power"` might be preferred for their frequentist operating characteristics, though though more research is needed on their general performance. `"multinom"` and `"poisson"` should only be used for comparison purposes or as an alternative interface to \pkg{boot}.
#'
#' @seealso
#' * [fwb.ci()] for calculating confidence intervals
#' * [summary.fwb()] for displaying output in a clean way
#' * [plot.fwb()] for plotting the bootstrap distributions
#' * [vcovFWB()] for estimating the covariance matrix of estimates using the FWB
#' * [set_fwb_wtype()] for an example of using weights other than the default exponential weights
#' * \pkgfun{boot}{boot} for the traditional bootstrap.
#'
#' See `vignette("fwb-rep")` for information on reproducibility.
#'
#' @references
#' Mammen, E. (1993). Bootstrap and Wild Bootstrap for High Dimensional Linear Models. *The Annals of Statistics*, 21(1). \doi{10.1214/aos/1176349025}
#'
#' Owen, A. B. (2025). Better bootstrap t confidence intervals for the mean. arXiv. \doi{10.48550/arXiv.2508.10083}
#'
#' Rubin, D. B. (1981). The Bayesian Bootstrap. *The Annals of Statistics*, 9(1), 130–134. \doi{10.1214/aos/1176345338}
#'
#' Xu, L., Gotwalt, C., Hong, Y., King, C. B., & Meeker, W. Q. (2020). Applications of the Fractional-Random-Weight Bootstrap. *The American Statistician*, 74(4), 345–358. \doi{10.1080/00031305.2020.1731599}
#'
#' @examplesIf rlang::is_installed("survival")
#' # Performing a Weibull analysis of the Bearing Cage
#' # failure data as done in Xu et al. (2020)
#' set.seed(123, "L'Ecuyer-CMRG")
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

#' @export
fwb <- function(data, statistic, R = 999, cluster = NULL, simple = NULL,
                wtype = getOption("fwb_wtype", "exp"), strata = NULL, drop0 = FALSE,
                verbose = TRUE, cl = NULL, ...) {

  bcall <- match.call()

  #Check data
  chk::chk_not_missing(data, "`data`")
  chk::chk_data(data)

  n <- nrow(data)
  if (is_null(n) || !chk::vld_count(n) || n < 1L) {
    .err("`data` must be present")
  }

  #Check statistic
  chk::chk_not_missing(statistic, "`statistic`")
  chk::chk_function(statistic)

  #Check R
  chk::chk_count(R)
  chk::chk_gt(R, 0)

  #Check cluster
  clus <- substitute(cluster)
  cluster <- eval(clus, data, parent.frame())

  if (is_not_null(cluster)) {
    .chk_atomic_vector(cluster)

    cluster <- factor(cluster)
    nc <- nlevels(cluster)
    cluster_numeric <- as.integer(cluster)
  }

  #Check strata
  strat <- substitute(strata)
  strata <- eval(strat, data, parent.frame())
  strata_to_use <- NULL

  if (is_not_null(strata)) {
    .chk_atomic_vector(strata)
    chk::chk_length(strata, n)

    strata_to_use <- factor(strata)

    if (is_not_null(cluster)) {
      cs <- unique(data.frame(cluster, strata_to_use))

      if (nrow(cs) != nlevels(cluster)) {
        .err("clusters must be completely nested within strata")
      }

      strata_to_use <- cs[[2L]]
    }
  }

  #Check verbose
  chk::chk_flag(verbose)

  #Check wtype
  gen_weights <- make_gen_weights(wtype)
  wtype <- .attr(gen_weights, "wtype")

  #Check drop0
  if (wtype %in% c("multinom", "poisson")) {
    chk::chk_scalar(drop0)
    chk::chk_logical(drop0)
  }
  else {
    drop0 <- FALSE
  }

  #Check simple
  if (is_not_null(simple)) {
    chk::chk_flag(simple)

    # if (simple && wtype == "multinom") {
    #   .err('`simple` cannot be `TRUE` when `wtype = "multinom"`')
    # }
  }
  else {
    simple <- wtype != "multinom"
  }

  #Process cl
  future.seed <- NULL
  if (is_not_null(cl)) {
    parallel_seed_set <- identical(RNGkind()[1L], "L'Ecuyer-CMRG")

    if (simple && !parallel_seed_set &&
        ((is.numeric(cl) && !isTRUE(all.equal(cl, 1))) || identical(cl, "future"))) {
      .wrn('`cl` was supplied but the random number generator is not suitable for parallelization. Set an appropriate seed using `set.seed(###, "L\'Ecuyer-CMRG")`, where ### is your favorite integer. See `?set.seed` for details')
    }

    if (identical(cl, "future")) {
      future.seed <- TRUE
    }
  }

  #Test fun
  test_w <- .set_class(rep.int(1, n), "fwb_internal_w")
  t0 <- try(call_statistic(statistic, data = data, ...,
                           .wi = test_w, drop0 = drop0))

  if (inherits(t0, "try-error")) {
    .err("There was an error running the function supplied to `statistic` on unit-weighted data. Error produced:\n\t",
         conditionMessage(.attr(t0, "condition")),
         tidy = FALSE)
  }

  if (!is.numeric(t0) || is_not_null(dim(t0))) {
    .err("the output of the function supplied to `statistic` must be a numeric vector")
  }

  if (anyNA(t0)) {
    .err("some estimates were returned as `NA` in the original sample")
  }

  random_statistic <- NULL
  if (simple) {
    t0_rep <- try(call_statistic(statistic, data = data, ...,
                                 .wi = test_w, drop0 = drop0))

    random_statistic <- !identical(t0_rep, t0)
  }

  if (is_null(names(t0))) {
    names(t0) <- paste0("t", seq_along(t0))
  }

  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    runif(1)
  }

  seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)

  w <- NULL
  if (is_null(cluster)) {
    if (simple) {
      FUN <- function(i) {
        wi <- drop(gen_weights(n, 1L, strata_to_use))
        call_statistic(statistic, data = data,
                       ...,
                       .wi = .set_class(wi, "fwb_internal_w"),
                       drop0 = drop0)
      }
    }
    else {
      w <- gen_weights(n, R, strata_to_use)
      FUN <- function(i) {
        call_statistic(statistic, data = data,
                       ...,
                       .wi = .set_class(w[i, ], "fwb_internal_w"),
                       drop0 = drop0)
      }
    }
  }
  else {
    if (simple) {
      FUN <- function(i) {
        cluster_w <- drop(gen_weights(nc, 1L, strata_to_use))
        wi <- cluster_w[cluster_numeric]
        call_statistic(statistic, data = data,
                       ...,
                       .wi = .set_class(wi, "fwb_internal_w"),
                       drop0 = drop0)
      }
    }
    else {
      cluster_w <- gen_weights(nc, R, strata_to_use)
      w <- cluster_w[, cluster_numeric, drop = FALSE]
      FUN <- function(i) {
        call_statistic(statistic, data = data,
                       ...,
                       .wi = .set_class(w[i, ], "fwb_internal_w"),
                       drop0 = drop0)
      }
    }
  }

  if (inherits(cl, "cluster")) {
    rlang::check_installed("parallel")

    parallel::clusterExport(cl, varlist = "call_statistic",
                            envir = asNamespace("fwb"))
  }

  opb <- pbapply::pboptions(type = if (verbose) "timer" else "none")
  on.exit(pbapply::pboptions(opb))

  #Run bootstrap
  t <- {
    if (is_not_null(future.seed))
      do.call("rbind", pbapply::pblapply(seq_len(R), FUN, cl = cl,
                                         future.seed = future.seed))
    else
      do.call("rbind", pbapply::pblapply(seq_len(R), FUN, cl = cl))
  }

  if (anyNA(t)) {
    .wrn("some estimates were returned as `NA`, which can cause problems in subsequent analyses")
  }

  colnames(t) <- names(t0)

  out <- list(t0 = t0,
              t = t,
              R = R,
              data = data,
              seed = seed,
              statistic = statistic,
              call = bcall,
              cluster = cluster,
              strata = strata,
              wtype = wtype)

  attr(out, "boot_type") <- "fwb"
  attr(out, "cl") <- cl
  attr(out, "simple") <- simple
  attr(out, "random_statistic") <- random_statistic

  class(out) <- c("fwb", "boot")

  out
}

#' @describeIn fwb
#' Print an `fwb` object
#'
#' @param x an `fwb` object; the output of a call to `fwb()`.
#' @param digits the number of significant digits to print
#' @param index the index or indices of the position of the quantity of interest in `x$t0` if more than one was specified in `fwb()`. Default is to print all quantities.
#'
#' @exportS3Method print fwb
print.fwb <- function(x, digits = getOption("digits", 3L), index = seq_len(ncol(x[["t"]])), ...) {
  index <- check_index(index, x[["t"]], several.ok = TRUE)
  t <- x[["t"]][, index, drop = FALSE]
  allNA <- apply(t, 2L, function(t_) all(is.na(t_)))
  ind1 <- index[allNA]
  index <- index[!allNA]
  t <- t[, !allNA, drop = FALSE]
  rn <- colnames(t)

  if (is_null(index)) {
    op <- NULL
  }
  else if (is_null(x$t0)) {
    op <- cbind(colMeans(t, na.rm = TRUE),
                apply(t, 2L, sd, na.rm = TRUE))
    dimnames(op) <- list(rn, c("mean", "std. error"))
  }
  else {
    op <- cbind(x$t0[index],
                colMeans(t, na.rm = TRUE) - x$t0[index],
                apply(t, 2L, sd, na.rm = TRUE))
    dimnames(op) <- list(rn, c("original", "bias", "std. error"))
  }

  if (is_null(x[["strata"]]) && is_null(x[["cluster"]])) {
    cat("FRACTIONAL WEIGHTED BOOTSTRAP\n")
  }
  else if (is_null(x[["strata"]]) && is_not_null(x[["cluster"]])) {
    cat("FRACTIONAL WEIGHTED CLUSTER BOOTSTRAP\n")
  }
  else if (is_not_null(x[["strata"]]) && is_null(x[["cluster"]])) {
    cat("STRATIFIED FRACTIONAL WEIGHTED BOOTSTRAP\n")
  }
  else {
    cat("STRATIFIED FRACTIONAL WEIGHTED CLUSTER BOOTSTRAP\n")
  }

  cat("\nCall:\n")
  dput(x$call, control = NULL)
  if (is_not_null(op)) {
    cat("\nBootstrap Statistics :\n")
    print(op, digits = digits)
  }

  if (is_not_null(ind1)) {
    for (j in ind1) {
      cat(sprintf("WARNING: All values of %s* are NA\n", colnames(x[["t"]])[j]))
    }
  }

  invisible(x)
}

check_statistic <- function(statistic) {

  statistic_args <- setdiff(names(formals(statistic)),
                            "...")

  if (length(statistic_args) < 2L) {
    .err("the function supplied to `statistic` must have at least two named arguments, the first corresponding to the dataset and the second corresponding to the weights")
  }

  forbidden_args <- setdiff(c(rlang::fn_fmls_names(fwb), rlang::fn_fmls_names(call_statistic)),
                            c("data", "..."))

  bad_args <- intersect(statistic_args, forbidden_args)

  if (is_not_null(bad_args)) {
    .err(sprintf("the function supplied to `statistic` cannot have arguments named %s",
                 word_list(bad_args, and.or = "or", quotes = TRUE)))
  }

  invisible(NULL)
}

make_gen_weights <- function(wtype) {
  chk::chk_string(wtype)
  wtype <- tolower(wtype)
  wtype <- match_arg(wtype, .w_types())

  fun <- switch(wtype,
                "exp" = function(n, R, strata = NULL) {
                  w <- matrix(rexp(n * R),
                              nrow = R, ncol = n, byrow = TRUE)

                  if (is_null(strata) || nlevels(strata) <= 1L) {
                    return(w / rowMeans(w))
                  }

                  for (s in levels(strata)) {
                    in_s <- which(strata == s)
                    w[, in_s] <- w[, in_s] / rowMeans(w[, in_s])
                  }

                  w
                },
                "poisson" = function(n, R, strata = NULL) {
                  matrix(rpois(n * R, 1),
                         nrow = R, ncol = n, byrow = TRUE)
                },
                "multinom" = function(n, R, strata = NULL) {
                  if (is_null(strata) || nlevels(strata) <= 1L) {
                    i <- sample.int(n, n * R, replace = TRUE)
                    dim(i) <- c(R, n)
                  }
                  else {
                    i <- matrix(NA_integer_, nrow = R, ncol = n)
                    for (s in levels(strata)) {
                      in_s <- which(strata == s)
                      n_s <- length(in_s)
                      i[, in_s] <- in_s[sample.int(n_s, n_s * R, replace = TRUE)]
                    }
                  }

                  t(apply(i, 1L, tabulate, n))
                },
                "mammen" = function(n, R, strata = NULL) {
                  sqrt5 <- sqrt(5)
                  w <- matrix((3 - sqrt5) / 2 + sqrt5 * rbinom(n * R, 1L, .5 - 1 / (2 * sqrt5)),
                              nrow = R, ncol = n, byrow = TRUE)

                  if (is_null(strata) || nlevels(strata) <= 1L) {
                    return(w / rowMeans(w))
                  }

                  for (s in levels(strata)) {
                    in_s <- which(strata == s)
                    w[, in_s] <- w[, in_s] / rowMeans(w[, in_s])
                  }

                  w
                },
                "beta" = function(n, R, strata = NULL) {
                  w <- matrix(4 * rbeta(n * R, .5, 1.5),
                              nrow = R, ncol = n, byrow = TRUE)

                  if (is_null(strata) || nlevels(strata) <= 1L) {
                    return(w / rowMeans(w))
                  }

                  for (s in levels(strata)) {
                    in_s <- which(strata == s)
                    w[, in_s] <- w[, in_s] / rowMeans(w[, in_s])
                  }

                  w
                },
                "power" = function(n, R, strata = NULL) {
                  w <- matrix((2 + sqrt(2)) * rbeta(n * R, sqrt(2) - 1, 1),
                              nrow = R, ncol = n, byrow = TRUE)

                  if (is_null(strata) || nlevels(strata) <= 1L) {
                    return(w / rowMeans(w))
                  }

                  for (s in levels(strata)) {
                    in_s <- which(strata == s)
                    w[, in_s] <- w[, in_s] / rowMeans(w[, in_s])
                  }

                  w
                })

  attr(fun, "wtype") <- wtype
  fun
}

call_statistic <- function(statistic, data, ..., .wi, drop0 = FALSE) {
  rlang::local_options(fwb_internal_w_env = rlang::current_env())

  if (!isFALSE(drop0)) {
    non0_wi <- which(.wi != 0)

    if (length(non0_wi) != length(.wi)) {
      if (isTRUE(drop0)) {
        .wi <- .wi[non0_wi]

        return(statistic(data[non0_wi, , drop = FALSE], .wi, ...))
      }

      if (is.na(drop0)) {
        is.na(.wi)[-non0_wi] <- TRUE
      }
    }
  }

  statistic(data, .wi, ...)
}

#' @exportS3Method stats::coef fwb
coef.fwb <- function(object, ...) {
  object[["t0"]]
}

#' @exportS3Method stats::vcov fwb
vcov.fwb <- function(object, ...) {
  stats::cov(object[["t"]])
}
