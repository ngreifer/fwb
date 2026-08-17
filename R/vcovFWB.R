#' Fractional Weighted Bootstrap Covariance Matrix Estimation
#'
#' `vcovFWB()` estimates the covariance matrix of model coefficient estimates using the fractional weighted bootstrap. It serves as a drop-in for [stats::vcov()] or \pkgfun{sandwich}{vcovBS}. Clustered covariances are can be requested.
#'
#' @inheritParams sandwich::vcovBS
#' @inheritParams fwb
#' @param x a fitted model object, such as the output of a call to `lm()` or `glm()`. The model object must result from a function that can be updated using [update()] and has a `weights` argument to input non-integer case weights.
#' @param R the number of bootstrap replications. When `wtype = "multinom"` and `R` is at least the number of distinct bootstrap samples of the clusters, every one of them is used exactly once instead and `R` is ignored; see [fwb()]. Default is 1000 (more is better but slower).
#' @param start `logical`; should `.coef(x)` be passed as `start` to the `update(x, weights = ...)` call? In case the model `x` is computed by some numeric iteration, this may speed up the bootstrapping. Default is `FALSE`.
#' @param wtype string; the type of weights to use. Allowable options include `"exp"` (the default), `"poisson"`, `"multinom"`, `"mammen"`, `"beta"`, and `"power"`. See [fwb()] for details. See [set_fwb_wtype()] to set a global default.
#' @param drop0 `logical`; when `wtype` is `"multinom"` or `"poisson"`, whether to drop units that are given weights of 0 from the model call in each iteration. If `TRUE`, the model will be called with an additional `subset` argument, filtering out units with weights of 0 (note this will overwrite any argument to `subset` in the original call). If `NA`, weights of 0 will be set to `NA` instead. Ignored for other `wtype`s because they don't produce 0 weights. Default is `FALSE`.
#' @param ... ignored.
#' @param fix `logical`; if `TRUE`, the covariance matrix is fixed to be positive semi-definite in case it is not.
#' @param use `character`; specification passed to [stats::cov()] for handling missing coefficients/parameters.
#' @param .coef a function used to extract the coefficients from each fitted model. Must return a numeric vector. By default, [`stats::coef`] is used, but `marginaleffects::get_coef` can be a more reliable choice for some models that have a non-standard `coef()` method, like that for `nnet::multinom()` models.
#'
#' @inherit sandwich::vcovBS return
#'
#' @details
#' `vcovFWB()` functions like other `vcov()`-like functions, such as those in the \pkg{sandwich} package, in particular, \pkgfun{sandwich}{vcovBS}, which implements the traditional bootstrap (and a few other bootstrap varieties for linear models). Sets of weights are generated as described in the documentation for [fwb()], and the supplied model is re-fit using those weights. When the fitted model already has weights, these are multiplied by the bootstrap weights.
#'
#' For `lm` objects, the model is re-fit using [.lm.fit()] for speed, and, similarly, `glm` objects are re-fit using [glm.fit()] (or whichever fitting method was originally used). For other objects, [update()] is used to populate the weights and re-fit the model (this assumes the fitting function accepts non-integer case weights through a `weights` argument). If a model accepts weights in some other way, [fwb()] should be used instead; `vcovFWB()` is inherently limited in its ability to handle all possible models. It is important that the original model was not fit using frequency weights (i.e., weights that allow one row of data to represent multiple full, identical, individual units) unless clustering is used.
#'
#' See \pkgfun{sandwich}{vcovBS} and \pkgfun{sandwich}{vcovCL} for more information on clustering covariance matrices, and see [fwb()] for more information on how clusters work with the fractional weighted bootstrap. When clusters are specified, each cluster is given a bootstrap weight, and all members of the cluster are given that weight; estimation then proceeds as normal. By default, when `cluster` is unspecified, each unit is considered its own cluster.
#'
#' ## Parallel Processing
#'
#' To speed up evaluation, parallel processing can be enabled. One way to do so is to supply an argument to `cl`. This can be either an integer (not available on Windows), a cluster object created by \pkgfun{parallel}{makeCluster}, or the string `"future"`. Another general way is to use functionality in the \CRANpkg{futurize} package, which is compatible with \pkg{fwb}. See `vignette("futurize-81-fwb", package = "futurize")` for details.
#'
#' Parallel processing does not change the results: the same seed gives the same covariance matrix whatever `cl` is set to, and calling [set.seed()] beforehand is all that is required to make a run reproducible. See `vignette("fwb-rep")` for details.
#'
#' When `cl` is a cluster object, \pkg{fwb} is attached on each worker so that the `w_*()` functions (see [w_mean()]) can be used in the model formula. Other packages the model needs must be attached by the user, e.g., with \pkgfun{parallel}{clusterEvalQ}.
#'
#' @seealso
#' * [fwb()] for performing the fractional weighted bootstrap on an arbitrary quantity
#' * [fwb.ci()] for computing nonparametric confidence intervals for `fwb` objects
#' * [summary.fwb()] for producing standard errors and confidence intervals for `fwb` objects
#' * \pkgfun{sandwich}{vcovBS} for computing covariance matrices using the traditional bootstrap (the fractional weighted bootstrap is also available but with limited options).
#'
#' @examplesIf rlang::is_installed("lmtest")
#' set.seed(123)
#' data("infert")
#' fit <- glm(case ~ spontaneous + induced, data = infert,
#'              family = "binomial")
#' lmtest::coeftest(fit, vcov. = vcovFWB, R = 200)
#'
#' @examplesIf rlang::is_installed("sandwich")
#' # Example from help("vcovBS", package = "sandwich")
#' data("PetersenCL", package = "sandwich")
#' m <- lm(y ~ x, data = PetersenCL)
#'
#' # Note: this is not to compare performance, just to
#' # demonstrate the syntax
#' cbind(
#'   "BS" = sqrt(diag(sandwich::vcovBS(m))),
#'   "FWB" = sqrt(diag(vcovFWB(m))),
#'   "BS-cluster" = sqrt(diag(sandwich::vcovBS(m, cluster = ~firm))),
#'   "FWB-cluster" = sqrt(diag(vcovFWB(m, cluster = ~firm)))
#' )
#'
#' # Using `wtype = "multinom"` exactly reproduces
#' # `sandwich::vcovBS()`
#' set.seed(11)
#' s <- sandwich::vcovBS(m, R = 200)
#' set.seed(11)
#' f <- vcovFWB(m, R = 200, wtype = "multinom")
#'
#' all.equal(s, f)
#' @examplesIf rlang::is_installed("nnet")
#' # Using a custom argument to `.coef`
#' set.seed(123)
#' data("infert")
#'
#' fit <- nnet::multinom(education ~ age, data = infert,
#'                       trace = FALSE)
#'
#' # vcovFWB(fit, R = 200) ## error
#' coef(fit) # coef() returns a matrix
#'
#' # Write a custom function to extract vector of
#' # coefficients (can also use marginaleffects::get_coef())
#' coef_multinom <- function(x) {
#'   p <- t(coef(x))
#'
#'   setNames(as.vector(p),
#'            paste(colnames(p)[col(p)],
#'                  rownames(p)[row(p)],
#'                  sep = ":"))
#' }
#' coef_multinom(fit) # returns a vector
#'
#' vcovFWB(fit, R = 200, .coef = coef_multinom)

#' @export
vcovFWB <- function(x, cluster = NULL, R = 1000, start = FALSE,
                    wtype = getOption("fwb_wtype", "exp"), drop0 = FALSE,
                    ..., fix = FALSE, use = "pairwise.complete.obs",
                    .coef = stats::coef,
                    verbose = FALSE, cl = NULL) {

  #Check arguments
  arg::arg_count(R)
  arg::arg_gt(R, 1)
  arg::arg_flag(start)
  arg::arg_flag(fix)
  arg::arg_string(use)
  arg::arg_flag(verbose)
  arg::arg_function(.coef)

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

  ## set up return value with correct dimension and names
  cf <- .coef(x)

  if (!is.numeric(cf) || length(dim(cf)) > 1L) {
    if (identical(.coef, eval(formals()[[".coef"]]))) {
      arg::err("the coefficients extracted using {.fun coef} from the supplied model are not in the form of a numeric vector; see the {.arg .coef} argument for {.fun vcovFWB}")
    }
    else {
      arg::err("the function supplied to {.arg .coef} must return a numeric vector")
    }
  }

  k <- length(cf)
  n <- nobs0(x)

  #One observation has exactly one possible bootstrap sample, so every replicate returns
  #the original coefficients and the covariance is all zeros. See `fwb()`.
  if (is_null(n) || n < 2L) {
    arg::err("the model must be fit to more than one unit to bootstrap; it was fit to {n %or% 0}")
  }

  rval <- matrix(0, nrow = k, ncol = k, dimnames = list(names(cf), names(cf)))

  if (is_null(cluster)) {
    cluster <- .attr(x, "cluster") %or% seq_len(n)
  }

  ## collect 'cluster' variables in a data frame
  if (inherits(cluster, "formula")) {
    cluster_tmp <- suppressWarnings(expand.model.frame(x, cluster, na.expand = FALSE))
    cluster <- model.frame(cluster, cluster_tmp, na.action = na.pass)
  }
  else {
    cluster <- as.data.frame(cluster)
  }

  ## handle omitted or excluded observations
  if (n != nrow(cluster) && is_not_null(x$na.action) &&
      inherits(x$na.action, c("exclude", "omit"))) {
    cluster <- cluster[-x$na.action, , drop = FALSE]
  }

  if (nrow(cluster) != n) {
    arg::err("number of observations in {.arg cluster} and {.fun nobs} do not match")
  }

  ## catch NAs in cluster -> need to be addressed in the model object by the user
  if (anyNA(cluster)) {
    arg::err("cannot handle {.val {NA}}s in {.arg cluster}: either refit the model without the {.val {NA}} observations in {.arg cluster} or impute the {.val {NA}}s")
  }

  ## every clustering variable must have something to resample
  for (j in seq_len(ncol(cluster))) {
    nc_j <- nlevels(factor(cluster[[j]]))

    if (nc_j < 2L) {
      arg::err("each variable in {.arg cluster} must have more than one cluster to bootstrap; {.val {names(cluster)[j]}} has {nc_j}")
    }
  }

  ## for multi-way clustering: set up interaction patterns
  p <- ncol(cluster)
  if (p > 1L) {
    clu <- lapply(seq_len(p), function(i) utils::combn(seq_len(p), i, simplify = FALSE))
    clu <- unlist(clu, recursive = FALSE)
    sgn <- (-1L)^(lengths(clu) + 1L)
    paste_ <- function(...) paste(..., sep = "_")
    for (i in (p + 1L):length(clu)) {
      cluster <- cbind(cluster, Reduce(paste_, unclass(cluster[, clu[[i]]]))) ## faster than: interaction()
    }
  }
  else {
    clu <- list(1)
    sgn <- 1
  }

  #Process cl. As in `fwb()`, each replicate's weights are drawn from a stream recorded
  #up front rather than from whatever the backend happens to hand the worker, so the
  #covariance does not depend on `cl`, on the worker count, or on the caller's RNG kind.
  future.seed <- NULL
  if (identical(cl, "future")) {
    future.seed <- TRUE
  }

  #A PSOCK worker has no packages attached beyond the defaults, so a model whose formula
  #uses `w_center()` and friends cannot be refit there without this. See `fwb()`.
  if (inherits(cl, "cluster")) {
    attach_on_workers(cl)
  }

  opb <- pbapply::pboptions(type = if (verbose) "timer" else "none")
  on.exit(pbapply::pboptions(opb))

  .env <- parent.frame()

  ## apply infrastructure for refitting models
  applyfun <- {
    if (is_not_null(future.seed))
      function(X, FUN, ...) pbapply::pblapply(X, FUN, ..., cl = cl, future.seed = future.seed)
    else
      function(X, FUN, ...) pbapply::pblapply(X, FUN, ..., cl = cl)
  }

  ## use starting values?
  start <- if (start && inherits(x, "glm")) .coef(x)

  ## bootstrap for each cluster dimension
  for (i in seq_along(clu)) {

    clu_i <- factor(cluster[[i]])

    ## exhaustive multinomial bootstrap, as in `fwb()`: once `R` reaches the number of
    ## distinct resamples, enumerate them instead of sampling. The count depends on the
    ## number of clusters, so it is decided per dimension -- one dimension can be
    ## exhaustive while a finer one is not, and each then uses its own number of
    ## replicates.
    R_i <- R
    w_all <- NULL

    if (wtype == "multinom") {
      n_resamples <- n_multinom_resamples(nlevels(clu_i))

      if (n_resamples <= R) {
        w_all <- gen_all_multinom_weights(nlevels(clu_i))
        R_i <- as.integer(n_resamples)
      }
    }

    ## one recorded stream per replicate; a fresh set per cluster dimension so the
    ## dimensions are drawn independently, as they are in `sandwich::vcovBS()`.
    ## Enumerated weights need no streams -- they consume no random numbers.
    seeds <- if (is_null(w_all)) make_stream_seeds(R_i)

    ## bootstrap fitting function
    bootfit <- make.bootfit(x, clu_i, start, drop0 = drop0,
                            gen_weights, .coef, .env, seeds, w_all)

    ## actually refit. `with_seed_preserved()` because `bootfit` installs a recorded
    ## stream into the global environment, which also sets the RNG kind
    cf <- with_seed_preserved(do.call("rbind", applyfun(seq_len(R_i), bootfit)))

    ## aggregate across cluster variables
    rval <- rval + sgn[i] * stats::cov(cf, use = use)
  }

  if (all_the_same(c(0, rval))) {
    arg::wrn("all variances and covariances are 0, indicating that the model failed to incorporate the bootstrapped weights")
    fix <- FALSE
  }

  ## check (and fix) if sandwich is not positive semi-definite
  if (fix) {
    eig <- eigen(rval, symmetric = TRUE)

    if (any(eig$values < 0)) {
      rval[] <- crossprod(sqrt(pmax(eig$values, 0)) * t(eig$vectors))
    }
  }

  rval
}

nobs0 <- function(x, ...) {
  if (inherits(x, "coxph")) {
    rval <- x[["n"]]
  }
  else {
    rval <- rlang::try_fetch(stats::nobs(x, ...),
                             error = function(e) NULL)
  }

  rval %or% NROW(residuals(x, ...))
}

make.bootfit <- function(fit, cli, start, drop0, gen_weights, .coef, .env, seeds = NULL,
                         w_all = NULL) {
  cli <- as.factor(cli)
  nc <- nlevels(cli)
  cluster_numeric <- as.integer(cli)

  #Replicate `j`'s cluster weights. `assign()` rather than package internals because this
  #closure is serialized to workers.
  #
  #`w_all` is the enumerated multinomial weight matrix, used when the exhaustive case
  #applies: replicate `j` is simply row `j`, with no random numbers involved. Otherwise
  #the weights are drawn from replicate `j`'s recorded stream, or -- with `seeds = NULL`,
  #which is how the tests build a `bootfit` directly -- from the ambient stream.
  gen_weights_j <- {
    if (is_not_null(w_all)) function(j) w_all[j, , drop = FALSE]
    else if (is.null(seeds)) function(j) gen_weights(nc, 1L)
    else function(j) {
      assign(".Random.seed", seeds[j, ], envir = globalenv())
      gen_weights(nc, 1L)
    }
  }
  special_coef <- !identical(.coef(fit), try(stats::coef(fit), silent = TRUE))
  bootfit <- NULL

  w0 <- weights(fit) %or% 1

  #The weights are computed at the length of the rows the fit actually used, but
  #`update()` re-evaluates the model call against the original data, so what it is handed
  #must be as long as that data. When `na.action` dropped rows, put each weight back at
  #the position it came from.
  #
  #Dropped rows are given a weight of 1. Any value would do: they are dropped again for
  #the same reason they were dropped the first time, so the value never reaches the
  #estimator. It does reach a `w_*()` call inside the model formula, which is evaluated
  #before `na.action` runs -- there, dropped rows count with weight 1. That matches how
  #base R treats `scale()` in a formula, which also sees every row.
  #
  #This is the inverse of the adjustment `vcovFWB()` makes to `cluster`, and rests on the
  #same assumption about what `na.action` indexes.
  na_act <- fit[["na.action"]]

  pad_weights <- {
    if (is_null(na_act) || !inherits(na_act, c("omit", "exclude"))) {
      identity
    }
    else {
      omitted <- unclass(na_act)

      function(w) {
        w_full <- rep.int(1, length(w) + length(omitted))
        w_full[-omitted] <- w
        w_full
      }
    }
  }

  if (!special_coef) {
    # Test if model uses w_*() functions by seeing if weighted
    # model.matrix differs from original; if so, need to use
    # general update instead of lm- or glm-specific.

    mm0 <- model.matrix(fit)

    if (is_not_null(mm0)) {
      #This refit always uses exponential weights, whatever `wtype` is, and never applies
      #`drop0`. All it needs is weights that are non-uniform and strictly positive.
      #
      #Non-uniform because that is the whole point: a `w_*()` term in the formula gives a
      #different model matrix under unequal weights than the unweighted fit did, and
      #comparing the two is how such a term is detected. Continuous and random rather
      #than a fixed pattern, because a regular one can coincide -- alternating weights of
      #1 and 2 reproduce the unweighted mean of any variable whose pattern has even
      #period, which would hide the very term being looked for.
      #
      #Strictly positive because some fitting functions reject a weight of 0 outright
      #(`survival::coxph()` does), so drawing from a `wtype` that produces zeros, or
      #applying `drop0`, would fail here on a call that was going to work. `drop0` is
      #irrelevant to the comparison anyway, and applying it changes the number of rows --
      #which the comparison would read as "the model matrix depends on the weights",
      #sending every model down the general path.
      probe_gen <- make_gen_weights("exp")

      probe_fit <- function() {
        rlang::local_options(fwb_internal_w_env = rlang::current_env())

        probe_w <- drop(probe_gen(nc, 1L))[cluster_numeric]

        .wi <- .set_class(pad_weights(probe_w), "fwb_internal_w")

        w <- .wi * pad_weights(w0)

        args <- list(fit,
                     weights = w,
                     evaluate = FALSE)

        if (is_not_null(start)) {
          args[["start"]] <- start
        }

        up <- do.call("update", args)

        utils::capture.output({
          suppressWarnings({
            up <- eval(up, envir = .env)
          })
        })

        model.matrix(up)
      }

      #A refit that fails tells us nothing, so fall back to the general path, which makes
      #no assumptions about the model. Better a slower covariance than none.
      #
      #`with_seed_preserved()` so that whether this probe runs, and how many random
      #numbers it draws, cannot shift the streams the replicates are drawn from.
      mm_test <- with_seed_preserved(
        rlang::try_fetch(probe_fit(),
                         error = function(e) NULL)
      )

      if (!isTRUE(all.equal(mm0, mm_test))) {
        special_coef <- TRUE
      }
    }
  }

  if (!special_coef && identical(class(fit), "lm")) {
    mf <- model.frame(fit)
    x <- model.matrix(fit)
    y <- model.response(mf) - (fit$offset %or% 0)

    w0 <- weights(fit) %or% 1

    p_x <- ncol(x)

    bootfit <- function(j, ...) {
      #Generate cluster weights, assign to units
      cluster_w <- drop(gen_weights_j(j))
      .wi <- cluster_w[cluster_numeric]

      .wi <- .wi * w0

      ws <- sqrt(.wi)

      z <- .lm.fit(x * ws, y = y * ws)

      cf <- z[["coefficients"]]

      #`.lm.fit()` is the raw decomposition: unlike `lm.fit()`, it leaves the coefficients
      #in *pivoted* order and leaves the aliased ones at whatever the decomposition
      #produced, which is usually 0. Both have to be undone for a replicate whose weights
      #leave the design rank-deficient, and the two do different damage.
      #
      #Without the `NA`s, an aliased coefficient enters `cov()` as a 0 that
      #`use = "pairwise.complete.obs"` cannot drop, so a replicate that estimated nothing
      #is averaged in as though it had.
      #
      #Without the reordering it is worse: `dqrdc2` moves negligible columns to the end,
      #so a coefficient is reported against the wrong parameter entirely. With three
      #columns and the second aliased, the pivot is `1 3 2` and the third column's
      #estimate lands in the second column's slot.
      #
      #The order here follows `lm.fit()`: mark the trailing (pivoted) positions `NA`
      #first, then restore the original column order.
      if (z[["rank"]] < p_x) {
        cf[seq.int(z[["rank"]] + 1L, p_x)] <- NA_real_
      }

      if (z[["pivoted"]]) {
        cf[z[["pivot"]]] <- cf
      }

      cf
    }
  }
  else if (!special_coef && identical(class(fit)[1L], "glm")) {
    x <- model.matrix(fit)
    y <- fit[["y"]]
    offset <- fit[["offset"]] %or% rep.int(0, nrow(x))

    w0 <- weights(fit) %or% 1

    if (is_null(fit[["method"]])) {
      fit.fun <- stats::glm.fit
    }
    else if (is.function(fit[["method"]])) {
      fit.fun <- fit[["method"]]
    }
    else if (is.character(fit[["method"]])) {
      fit.fun <- get0(fit[["method"]], envir = environment(fit[["terms"]]),
                      mode = "function")
      if (is_null(fit.fun)) {
        arg::err("the {.arg method} used to fit the original model ({.val {fit[['method']]}}) is unavailable")
      }
    }
    else {
      arg::err("unrecognized fitting method; the model cannot be re-fit")
    }

    bootfit <- function(j, ...) {
      #Generate cluster weights, assign to units
      cluster_w <- drop(gen_weights_j(j))
      .wi <- cluster_w[cluster_numeric]

      .wi <- .wi * w0

      #`drop0 = NA` means "set weights of 0 to `NA`", which everywhere else works because
      #the model's own `na.action` then drops those rows. There is no `na.action` in a
      #direct call to the fitting function, and it cannot take `NA` weights at all, so
      #here the two settings coincide: dropping the rows *is* what setting their weights
      #to `NA` accomplishes. The estimates are identical either way.
      if (!isFALSE(drop0)) {
        zero_w <- .wi == 0

        if (any(zero_w)) {
          .wi <- .wi[!zero_w]

          x <- x[!zero_w, , drop = FALSE]
          y <- y[!zero_w]
          offset <- offset[!zero_w]
        }
      }

      safe.glm.fit(fit.fun, x = x, y = y,
                   weights = .wi,
                   start = start,
                   offset = offset,
                   family = fit$family,
                   control = fit$control,
                   intercept = .attr(fit$terms, "intercept") > 0)$coefficients
    }
  }
  else {
    w0 <- weights(fit) %or% 1

    bootfit <- function(j, ...) {
      rlang::local_options(fwb_internal_w_env = rlang::current_env())

      #Generate cluster weights, assign to units
      cluster_w <- drop(gen_weights_j(j))
      .wi <- cluster_w[cluster_numeric]

      .wi <- .set_class(pad_weights(.wi), "fwb_internal_w")

      w <- .wi * pad_weights(w0)

      args <- list(fit,
                   weights = w,
                   evaluate = FALSE)

      if (isTRUE(drop0)) {
        args[["subset"]] <- w > 0
      }
      else if (is.na(drop0)) {
        is.na(args[["weights"]][w == 0]) <- TRUE
      }

      if (is_not_null(start)) {
        args[["start"]] <- start
      }

      #Swallows chatty fitting functions. Costs ~18us per replicate, which is 1-3%
      #of a `coxph` refit but up to 8% of a cheap one (an `lm()` pulled onto this
      #path by a `w_*()` term). A single `sink()` around the whole `applyfun()`
      #would recover that, but it would also swallow the progress bar and leak an
      #open sink if a replicate errored; see dev/review-2026-08.md.
      utils::capture.output({
        up <- do.call("update", args)

        up <- eval(up, envir = .env)
      })

      .coef(up)
    }
  }

  bootfit
}

safe.glm.fit <- function(fit.fun, ...) {
  withCallingHandlers({
    fit.fun(...)
  },
  warning = function(w) {
    if (conditionMessage(w) != "non-integer #successes in a binomial glm!" &&
        !startsWith(conditionMessage(w), "non-integer x =")) {
      arg::wrn("{w}", immediate = FALSE)
    }
    invokeRestart("muffleWarning")
  })
}
