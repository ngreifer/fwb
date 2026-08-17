#' Fractional Weighted Bootstrap
#'
#' `fwb()` implements the fractional (random) weighted bootstrap, also known as the Bayesian bootstrap. Rather than resampling units to include in bootstrap samples, weights are drawn to be applied to a weighted estimator.
#'
#' @param data the dataset used to compute the statistic; must contain more than one unit
#' @param statistic a function, which, when applied to `data`, returns a vector containing the statistic(s) of interest. The function should take at least two arguments; the first argument should correspond to the dataset and the second argument should correspond to a vector of weights. Any further arguments can be passed to `statistic` through the `...` argument, but they cannot share a name with an argument of `fwb()`, which would claim them instead. These requirements are checked unless `statistic` has a `...` argument, in which case no assumption is made about what it accepts.
#' @param R the number of bootstrap replicates. Default is 999 but more is always better. For the percentile bootstrap confidence interval to be exact, it can be beneficial to use one less than a multiple of 100. When `wtype = "multinom"` and `R` is at least the number of distinct bootstrap samples of `data`, it is ignored; see the description of `"multinom"` below.
#' @param cluster optional; a vector containing cluster membership. If supplied, will run the cluster bootstrap. Must contain more than one cluster. See Details. Evaluated first in `data` and then in the global environment.
#' @param simple `logical`; if `TRUE`, weights will be generated on-the-fly in each bootstrap replication; if `FALSE`, all weights will be generated at once and then supplied to `statistic`. The default (`NULL`) sets to `FALSE` if `wtype = "multinom"` and to `TRUE` otherwise.
#' @param wtype string; the type of weights to use. Allowable options include `"exp"` (the default), `"poisson"`, `"multinom"`, `"mammen"`, `"beta"`, and `"power"`. See Details. See [set_fwb_wtype()] to set a global default.
#' @param strata optional; a vector containing stratum membership for stratified bootstrapping. If supplied, will essentially perform a separate bootstrap within each level of `strata`. This does not affect results when `wtype = "poisson"`.
#' @param drop0 `logical`; when `wtype` is `"multinom"` or `"poisson"`, whether to drop units that are given weights of 0 from the dataset and weights supplied to `statistic` in each iteration. If `NA`, weights of 0 will be set to `NA` instead. Ignored for other `wtype`s because they don't produce 0 weights. Default is `FALSE`.
#' @param verbose `logical`; whether to display a progress bar. The default value, `NULL`, is `FALSE` when parallelization is used (see `cl` below) and `TRUE` otherwise. Alternatively, the \CRANpkg{progressify} package can be used to incorporate a \pkg{progressr} progress bar.
#' @param cl a cluster object created by \pkgfun{parallel}{makeCluster}, an integer to indicate the number of child-processes (integer values are ignored on Windows) for parallel evaluations, or the string `"future"` to use a `future` backend. See the `cl` argument of \pkgfun{pbapply}{pblapply} for details. If `NULL`, no parallelization will take place. Alternatively, the \CRANpkg{futurize} package can be used to incorporate a `future` backend. See the section "Parallel Processing" in Details.
#' @param ... other arguments passed to `statistic`.
#'
#' @returns
#' An `<fwb>` object, which also inherits from `<boot>`, with the following components:
#'
#' \item{t0}{The observed value of `statistic` applied to `data` with uniform weights.}
#' \item{t}{A matrix with `R` rows, each of which is a bootstrap replicate of the result of calling `statistic`.}
#' \item{R}{The value of `R` as passed to `fwb()`, except in the exhaustive multinomial case described under `"multinom"` below, where it is the number of distinct bootstrap samples.}
#' \item{data}{The `data` as passed to `fwb()`.}
#' \item{statistic}{The function `statistic` as passed to `fwb()`.}
#' \item{call}{The original call to `fwb()`.}
#' \item{cluster}{The vector passed to `cluster`, if any.}
#' \item{strata}{The vector passed to `strata`, if any.}
#' \item{wtype}{The type of weights used as determined by the `wtype` argument.}
#'
#' The seed (when `simple = FALSE`) or seeds (when `simple = TRUE`) used to re-generate the weights is stored in the `"seeds"` attribute of the returned object. In the exhaustive multinomial case described under `"multinom"` below, no random numbers are used and the `"exhaustive"` attribute is `TRUE`.
#'
#' `<fwb>` objects have [coef()] and [vcov()] methods, which extract the `t0` component and covariance of the `t` components, respectively.
#'
#' @details
#' `fwb()` implements the fractional weighted bootstrap and is meant to function as a drop-in for `boot::boot(., stype = "f")` (i.e., the usual bootstrap but with frequency weights representing the number of times each unit is drawn). In each bootstrap replication, when `wtype = "exp"` (the default), the weights are sampled from independent exponential distributions with rate parameter 1 and then normalized to have a mean of 1, equivalent to drawing the weights from a Dirichlet distribution. Other weights are allowed as determined by the `wtype` argument (see below for details). The function supplied to `statistic` must incorporate the weights to compute a weighted statistic. For example, if the output is a regression coefficient, the weights supplied to the `w` argument of `statistic` should be supplied to the `weights` argument of `lm()`. These weights should be used any time frequency weights would be, since they are meant to function like frequency weights (which, in the case of the traditional bootstrap, would be integers). Unfortunately, there is no way for `fwb()` to know whether you are using the weights correctly, so care should be taken to ensure weights are correctly incorporated into the estimator.
#'
#' When fitting binomial regression models (e.g., logistic) using [glm()], it may be useful to change the `family` to a "quasi" variety (e.g., [quasibinomial()]) to avoid a spurious warning about "non-integer #successes".
#'
#' The cluster/block bootstrap can be requested by supplying a vector of cluster membership to `cluster`. Rather than generating a weight for each unit, a weight is generated for each cluster and then applied to all units in that cluster.
#'
#' Bootstrapping can be performed within strata by supplying a vector of stratum membership to `strata`. This essentially rescales the weights within each stratum to have a mean of 1, ensuring that the sum of weights in each stratum is equal to the stratum size. For multinomial weights, using strata is equivalent to drawing samples with replacement from each stratum. Strata do not affect bootstrapping when using Poisson weights.
#'
#' The bootstrap weights depend only on the state of the random number generator when `fwb()` is called, so calling [set.seed()] beforehand is all that is required to make a run reproducible. In particular, the weights do not depend on `cl`, on the number of workers, on `verbose`, or on whether `statistic` itself draws random numbers, and no particular [RNGkind()] is required. See `vignette("fwb-rep")` for details.
#'
#' Note that `simple` does affect the weights: `simple = FALSE` draws them in a single call before `statistic` is applied, whereas `simple = TRUE` draws each replicate's weights from a random number stream reserved for that replicate. Both are reproducible, but they are not interchangeable, so you will get different results depending on whether `simple` is `TRUE` or `FALSE`. It's recommended to leave `simple` at its default value.
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
#'
#' Unlike the other weight types, the multinomial weights have only finitely many possible values: there are \eqn{n^n} distinct bootstrap samples of \eqn{n} units (or \eqn{\prod_s n_s^{n_s}} with strata, and computed at the cluster level when `cluster` is supplied). When `R` is at least that many, `fwb()` enumerates every one of them exactly once instead of sampling, sets `R` to that number, and reports this in a message. The resulting estimates are then the exact bootstrap distribution rather than a sample from it, and they no longer depend on the seed. This is the one respect in which `wtype = "multinom"` differs from `boot::boot(., stype = "f")`, which always samples; it only applies to very small samples, as \eqn{n = 6} already has 46656 distinct samples.
#' }
#' \item{`"poisson"`}{
#' Draws integer weights from a Poisson distribution with \eqn{\lambda = 1} using [rpois()]. This is an alternative to the multinomial weights that yields similar estimates (especially as the sample size grows) but can be faster. Note `strata` is ignored when using `"poisson"`.
#' }
#' \item{`"mammen"`}{
#' Draws weights from a modification of the distribution described by Mammen (1983) for use in the wild bootstrap. These positive weights have a mean, variance, and skewness of 1, making them second-order accurate (in contrast to the usual exponential weights, which are only first-order accurate). The weights \eqn{w} are drawn such that \eqn{P(w=(3+\sqrt{5})/2)=(\sqrt{5}-1)/2\sqrt{5}} and \eqn{P(w=(3-\sqrt{5})/2)=(\sqrt{5}+1)/2\sqrt{5}} as described by Owen (2025). The weights are scaled to have a mean of 1 within each stratum (or in the full sample if `strata` is not supplied).
#' }
#' \item{`"beta"`}{
#' Draws weights from a \eqn{\text{Beta}(1/2, 3/2)} distribution using [rbeta()] as described by Owen (2025). These positive weights have a mean, variance, and skewness of 1 when scaled by a factor of 4, making them second-order accurate. The weights are scaled to have a mean of 1 within each stratum (or in the full sample if `strata` is not supplied).
#' }
#' \item{`"power"`}{
#' Draws weights from a \eqn{\text{Beta}(\sqrt{2} - 1, 1)} distribution using [rbeta()] as described by Owen (2025). These positive weights have a mean and variance of 1 and skewness of \eqn{2(\sqrt{2} - 1)} when scaled by a factor of \eqn{2+\sqrt{2}}. The weights are scaled to have a mean of 1 within each stratum (or in the full sample if `strata` is not supplied).
#' }
#' }
#'
#' `"exp"` is the default due to it being the formulation described in Xu et al. (2020) and in the most formulations of the Bayesian bootstrap; it should be used if one wants to remain in line with these guidelines or to maintain a Bayesian flavor to the analysis, whereas other distributions might be preferred for their frequentist operating characteristics, though more research is needed on their general performance. Owen (2025) recommends `"beta"` and `"power"`, as these provided close to nominal confidence interval coverage without excessively large intervals in the context of estimating the mean in a small sample. `"multinom"` and `"poisson"` should primarily be used for comparison purposes or as an alternative interface to \pkg{boot}.
#'
#' ## Parallel Processing
#'
#' To speed up evaluation, parallel processing can be enabled. One way to do so is to supply an argument to `cl`. This can be either an integer (not available on Windows), a `<cluster>` object created by \pkgfun{parallel}{makeCluster}, or the string `"future"` (combined with setting a parallelization scheme with \pkgfun{future}{plan}. Another general way is to use functionality in the \CRANpkg{futurize} package, which is compatible with \pkg{fwb}. See `vignette("futurize-81-fwb", package = "futurize")` for details.
#'
#' Parallel processing does not change the results: the same seed gives the same bootstrap replicates whatever `cl` is set to. See `vignette("fwb-rep")` for details.
#'
#' When `cl` is a `<cluster>` object, \pkg{fwb} is attached on each worker so that the `w_*()` functions (e.g., [w_mean()], [w_std()]) can be used inside `statistic`. Other packages needed by `statistic` must be attached by the user, e.g., with \pkgfun{parallel}{clusterEvalQ}.
#'
#' When parallel processing is used, a \CRANpkg{progressr} progress bar via the \CRANpkg{progressify} package requires a \pkg{future} backend (i.e., `cl = "future"` or the `futurize()` pipe); with other values of `cl` no progress is reported. The \pkg{pbapply} progress bar controlled by `verbose` works with all of them, but may slow down evaluation and so is not recommended. \pkg{progressify} progress bars work as normal with no parallelization.
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
#' set.seed(123)
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
                verbose = NULL, cl = NULL, ...) {

  bcall <- match.call()

  #Check data
  arg::arg_supplied(data)
  arg::arg_data.frame(data)

  n <- nrow(data)
  if (is_null(n) || !rlang::is_scalar_integerish(n) || n < 1L) {
    arg::err("{.arg data} must be present")
  }

  #A single unit has exactly one possible bootstrap sample -- itself -- so every
  #replicate is the original estimate and every standard error is 0. Nothing in the
  #package means anything in that case, so refuse rather than return zeros.
  if (n < 2L) {
    arg::err("{.arg data} must have more than one unit to bootstrap; it has {n}")
  }

  #Check statistic
  arg::arg_supplied(statistic)
  arg::arg_function(statistic)
  check_statistic(statistic)
  .statistic <- substitute(statistic)

  #Check R
  arg::arg_count(R)
  arg::arg_gt(R, 0)

  #Check cluster
  clus <- substitute(cluster)
  cluster <- eval(clus, data, parent.frame())

  arg::when_not_null(cluster,
                     arg::arg_vector,
                     arg::arg_length(n))

  if (is_not_null(cluster)) {
    cluster <- factor(cluster)
    nc <- nlevels(cluster)
    cluster_numeric <- as.integer(cluster)

    #The cluster bootstrap resamples clusters, so one cluster is the `n < 2` case above:
    #the only possible sample is the original one.
    if (nc < 2L) {
      arg::err("{.arg cluster} must have more than one cluster to bootstrap; it has {nc}")
    }
  }

  #Check strata
  .strata <- substitute(strata)
  strata <- eval(.strata, data, parent.frame())
  strata_to_use <- NULL

  arg::when_not_null(strata,
                     arg::arg_vector,
                     arg::arg_length(n))

  if (is_not_null(strata)) {
    strata_to_use <- factor(strata)

    if (is_not_null(cluster)) {
      cs <- unique(data.frame(cluster, strata_to_use))

      if (nrow(cs) != nlevels(cluster)) {
        arg::err("clusters must be completely nested within strata")
      }

      strata_to_use <- cs[[2L]]
    }
  }

  #Check verbose
  arg::when_not_null(verbose, arg::arg_flag)

  if (is_null(verbose)) {
    verbose <- (guess_num_workers(cl) == 1)
  }

  #Check wtype
  gen_weights <- make_gen_weights(wtype)
  wtype <- .attr(gen_weights, "wtype")

  #Check drop0
  if (wtype %in% c("multinom", "poisson")) {
    arg::arg_or(drop0,
                arg::arg_is_NA,
                arg::arg_flag)
  }
  else {
    drop0 <- FALSE
  }

  #Check simple
  arg::when_not_null(simple, arg::arg_flag)

  if (is_null(simple)) {
    simple <- wtype != "multinom"
  }

  #Weights are drawn at the cluster scale when there are clusters, then expanded onto
  #units; `expand` is empty when there are none.
  n_w <- if (is_null(cluster)) n else nc
  expand <- if (is_null(cluster)) integer(0L) else cluster_numeric

  #Exhaustive multinomial bootstrap.
  #
  #The multinomial weights have only finitely many possible values, so once `R` reaches
  #that many, sampling at random is strictly worse than enumeration: it revisits some
  #resamples and misses others, and the resulting distribution is a sample from the
  #bootstrap distribution rather than the bootstrap distribution itself. Enumerating
  #instead makes the estimates exact and makes them deterministic -- the same for every
  #seed -- so `R` is reduced to the number of distinct resamples and no more are drawn.
  #
  #`simple` becomes meaningless here and is set to `FALSE` accordingly: there is no
  #randomness left to reserve a stream for, and the whole matrix is built at once.
  #
  #This is the one place `wtype = "multinom"` departs from `boot::boot(., stype = "f")`,
  #which always samples. It reaches only very small samples: `n = 6` already has 46656
  #resamples, so with the default `R` nothing above `n = 4` is affected.
  exhaustive <- FALSE

  if (wtype == "multinom") {
    n_resamples <- n_multinom_resamples(n_w, strata_to_use)

    #Reached when every stratum holds a single unit: the only possible resample of a
    #stratum of one is itself, so all weights are 1 in every replicate.
    if (n_resamples <= 1) {
      arg::err(c("there is only one possible multinomial bootstrap sample of these data, so there is nothing to bootstrap.",
                 i = "this happens when every level of {.arg strata} contains a single unit"))
    }

    if (n_resamples <= R) {
      exhaustive <- TRUE
      simple <- FALSE

      if (n_resamples < R) {
        arg::msg("there are only {n_resamples} distinct multinomial bootstrap samples of these data; using all of them and ignoring {.arg R}")
      }

      R <- as.integer(n_resamples)
    }
  }

  #Process cl. When `simple = TRUE` the weights are seeded per replicate (see below),
  #so the backend never supplies the randomness for them and no particular RNG kind is
  #required of the caller. `simple = FALSE` draws the whole matrix here in the parent,
  #so it needs nothing from the backend either.
  #
  #`future.seed = TRUE` is still passed so that `future.apply` does not warn about a
  #future touching the RNG; whatever seed it installs is overwritten before any draw.
  future.seed <- NULL
  if (identical(cl, "future")) {
    future.seed <- TRUE
  }

  #Record one stream per replicate before `statistic` has been called even once, so that
  #the weights depend on the seed, `R`, and the RNG kind and on nothing else. Deriving
  #them after the unit-weight call below would let a `statistic` that draws random
  #numbers shift every subsequent weight.
  seeds <- {
    if (simple) make_stream_seeds(R)
    else NULL
  }

  #Test fun
  test_w <- .set_class(rep.int(1, n), "fwb_internal_w")
  t0 <- rlang::try_fetch(
    call_statistic(statistic, data = data, ...,
                   .wi = test_w, drop0 = drop0),
    error = function(e) {
      if (rlang::is_call(e$call, "statistic")) {
        e$call[[1L]] <- .statistic
      }

      arg::err("there was a problem running the function supplied to {.arg statistic} on unit-weighted data",
               parent = e)
    }
  )

  if (!is.numeric(t0) || is_not_null(dim(t0))) {
    arg::err("the output of the function supplied to {.arg statistic} must be a numeric vector")
  }

  if (anyNA(t0)) {
    arg::err("some estimates were returned as {.val {NA}} in the original sample")
  }

  if (is_null(names(t0))) {
    names(t0) <- paste0("t", seq_along(t0))
  }

  #Two ways of drawing the weights, and each records what `fwb.array()` needs to
  #recover them without re-running `statistic`.
  #
  #`simple = FALSE` draws the whole matrix here, so `seed` -- the state just before
  #that one call -- is enough. This also keeps `wtype = "multinom"` byte-identical to
  #`boot::boot(., stype = "f")`, which consumes its stream the same way.
  #
  #`simple = TRUE` draws inside each replicate, where the backend decides the order and
  #whether the stream is split. Rather than trying to replay that, each replicate gets
  #its own recorded stream: the weights become a function of the replicate index alone,
  #identical for every `cl` and `verbose`, and recoverable exactly even when
  #`statistic` itself draws random numbers.
  #
  #The exhaustive multinomial case takes the `simple = FALSE` branch by construction and
  #consumes no random numbers at all, so the recorded seed is irrelevant to it.

  w <- NULL
  if (simple) {
    #`assign()` rather than the `set_stream_seed()` helper, and `length(expand)` rather
    #than `is_not_null()`: this closure is serialized to workers, so the fewer package
    #internals its body has to resolve there, the fewer ways it can fail.
    FUN <- function(i) {
      assign(".Random.seed", seeds[i, ], envir = globalenv())

      wi <- drop(gen_weights(n_w, 1L, strata_to_use))

      if (length(expand) > 0L) {
        wi <- wi[expand]
      }

      call_statistic(statistic, data = data, ...,
                     .wi = .set_class(wi, "fwb_internal_w"),
                     drop0 = drop0)
    }
  }
  else {
    seeds <- get(".Random.seed", envir = globalenv(), inherits = FALSE)

    w <- {
      if (exhaustive) gen_all_multinom_weights(n_w, strata_to_use)
      else gen_weights(n_w, R, strata_to_use)
    }

    if (length(expand) > 0L) {
      w <- w[, expand, drop = FALSE]
    }

    FUN <- function(i) {
      call_statistic(statistic, data = data, ...,
                     .wi = .set_class(w[i, ], "fwb_internal_w"),
                     drop0 = drop0)
    }
  }

  #A PSOCK worker starts with no packages attached beyond the defaults, and
  #`statistic` resolves names through its own enclosing environment -- usually the
  #global environment, which is empty on the worker. So a `statistic` calling
  #`w_mean()` and friends fails with "could not find function". Attaching the package
  #on each worker is what makes the `w_*()` helpers usable with a `cluster`, as they
  #already are with forking and with a `future` backend.
  if (inherits(cl, "cluster")) {
    attach_on_workers(cl)
  }

  opb <- pbapply::pboptions(type = if (verbose) "timer" else "none")
  on.exit(pbapply::pboptions(opb))

  #Run bootstrap.
  #
  #`with_seed_preserved()` matters when `simple = TRUE` and evaluation is sequential:
  #`FUN` installs a recorded stream into the global environment, which also switches the
  #RNG kind. Restoring afterwards leaves the caller's kind and stream where
  #`make_stream_seeds()` left them -- one draw on from where they started, so successive
  #calls to `fwb()` use different streams.
  run <- function() {
    if (is_not_null(future.seed))
      do.call("rbind", pbapply::pblapply(seq_len(R), FUN, cl = cl,
                                         future.seed = future.seed))
    else
      do.call("rbind", pbapply::pblapply(seq_len(R), FUN, cl = cl))
  }

  t <- {
    if (simple) with_seed_preserved(run())
    else run()
  }

  if (anyNA(t)) {
    arg::wrn("some estimates were returned as {.val {NA}}, which can cause problems in subsequent analyses")
  }

  colnames(t) <- names(t0)

  out <- list(t0 = t0,
              t = t,
              R = R,
              data = data,
              statistic = statistic,
              call = bcall,
              cluster = cluster,
              strata = strata,
              wtype = wtype)

  attr(out, "boot_type") <- "fwb"
  attr(out, "simple") <- simple
  attr(out, "seeds") <- seeds
  attr(out, "exhaustive") <- exhaustive

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

#Validate `statistic`, without assuming it was written by the user.
#
#A function with `...` is exempt from every check below. That is not a loophole, it is
#the point: `fwb()` calls `statistic(data, w, ...)`, so a function with `...` can accept
#whatever it is handed, and nothing about its formals can be shown to be wrong. It is
#also what makes this safe for other packages to build on. Anything that wraps
#`statistic` -- *progressify*'s transpiler, or something written later -- has to forward
#arguments to the function it wraps, and therefore has `...`. Keying the exemption on
#`...` rather than on a list of known packages means such wrappers pass without `fwb`
#having to know they exist.
check_statistic <- function(statistic) {

  statistic_args <- rlang::fn_fmls_names(statistic) %or% character(0L)

  if ("..." %in% statistic_args) {
    return(invisible(NULL))
  }

  #`fwb()` supplies the data and the weights positionally, so a function that cannot
  #take two arguments can never receive the weights at all.
  if (length(statistic_args) < 2L) {
    arg::err(c("the function supplied to {.arg statistic} must accept at least two arguments, the first the dataset and the second the weights.",
               i = "it has {length(statistic_args)}: {.or {.arg {statistic_args}}}"))
  }

  #Further arguments reach `statistic` only through `fwb()`'s own `...`, and are matched
  #by name, so an argument sharing a name with one of `fwb()`'s can never be given a
  #value -- `fwb()` claims it first. The first two are positional, so their names are
  #unconstrained.
  forbidden_args <- setdiff(rlang::fn_fmls_names(fwb), c("data", "..."))

  bad_args <- intersect(statistic_args[-(1:2)], forbidden_args)

  if (is_not_null(bad_args)) {
    arg::err(c("the function supplied to {.arg statistic} cannot have argument{?s} named {.or {.arg {bad_args}}}.",
               i = "{cli::qty(length(bad_args))}{.fun fwb} would claim {?it/them} as its own, so {?it/they} could never be supplied through {.arg ...}"))
  }

  invisible(NULL)
}

#Number of distinct multinomial bootstrap samples.
#
#The multinomial bootstrap draws `m` units with replacement from `m`, so it has exactly
#`m^m` possible outcomes. Strata are drawn independently, so the count is the product
#over strata. Returned as a double because the count leaves integer range almost at once
#(`m = 10` is already 1e10) and every caller only compares it against `R`.
n_multinom_resamples <- function(n, strata = NULL) {
  sizes <- {
    if (is_null(strata) || nlevels(strata) <= 1L) n
    else tabulate(strata, nlevels(strata))
  }

  sizes <- sizes[sizes > 0L]

  prod(as.double(sizes)^sizes)
}

#Every ordered resample of `seq_len(m)`, one per row of an `m^m x m` matrix: the rows are
#the length-`m` numbers in base `m`, so position `p` cycles with period `m^(p - 1)`.
all_resamples <- function(m) {
  N <- m^m

  i <- matrix(0L, nrow = N, ncol = m)

  for (p in seq_len(m)) {
    i[, p] <- rep(seq_len(m), each = m^(p - 1L), length.out = N)
  }

  i
}

#Every multinomial weight matrix that could be drawn, in a fixed order.
#
#Deterministic -- no random numbers are drawn -- which is what lets `fwb.array()` and
#`vcovFWB()` reproduce the weights without recording anything.
#
#The enumeration is over ordered resamples rather than over distinct weight vectors, so
#each weight vector appears as often as its multinomial probability requires: at `n = 2`
#the four resamples give weights `(2,0)`, `(1,1)`, `(1,1)`, `(0,2)`, which is the exact
#multinomial distribution. Enumerating the three distinct vectors once each would not be.
gen_all_multinom_weights <- function(n, strata = NULL) {
  if (is_null(strata) || nlevels(strata) <= 1L) {
    return(tabulate_rows(all_resamples(n), n))
  }

  in_s <- lapply(levels(strata), function(s) which(strata == s))
  in_s <- in_s[lengths(in_s) > 0L]

  w_s <- lapply(in_s, function(u) tabulate_rows(all_resamples(length(u)), length(u)))

  n_s <- vapply(w_s, nrow, integer(1L))

  N <- prod(n_s)

  w <- matrix(0L, nrow = N, ncol = n)

  #Mixed radix over the strata: stratum `k` varies once every `before` rows, so every
  #combination of per-stratum resamples appears exactly once.
  before <- 1

  for (k in seq_along(in_s)) {
    rows <- ((seq_len(N) - 1L) %/% before) %% n_s[k] + 1L

    w[, in_s[[k]]] <- w_s[[k]][rows, , drop = FALSE]

    before <- before * n_s[k]
  }

  w
}

make_gen_weights <- function(wtype) {
  wtype <- arg::match_arg(wtype, .w_types())

  fun <- switch(wtype,
                "exp" = function(n, R, strata = NULL) {
                  w <- matrix(rexp(n * R),
                              nrow = R, ncol = n, byrow = TRUE)

                  if (is_null(strata) || nlevels(strata) <= 1L) {
                    return(w / rowMeans(w))
                  }

                  for (s in levels(strata)) {
                    in_s <- which(strata == s)
                    w[, in_s] <- w[, in_s, drop = FALSE] / rowMeans(w[, in_s, drop = FALSE])
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

                  tabulate_rows(i, n)
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
                    w[, in_s] <- w[, in_s, drop = FALSE] / rowMeans(w[, in_s, drop = FALSE])
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
                    w[, in_s] <- w[, in_s, drop = FALSE] / rowMeans(w[, in_s, drop = FALSE])
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
                    w[, in_s] <- w[, in_s, drop = FALSE] / rowMeans(w[, in_s, drop = FALSE])
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
        is.na(.wi[-non0_wi]) <- TRUE
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
