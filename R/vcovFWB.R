#' Fractional Weighted Bootstrap Covariance Matrix Estimation
#'
#' `vcovFWB()` estimates the covariance matrix of model coefficient estimates using the fractional weighted bootstrap. It serves as a drop-in for `stats::vcov()` or `sandwich::vcovBS()`. Clustered covariances are can be requested.
#'
#' @inheritParams sandwich::vcovBS
#' @inheritParams fwb
#' @param x a fitted model object, such as the output of a call to `lm()` or `glm()`. The model object must result from a function that can be updated using [update()] and has a `weights` argument to input non-integer case weights.
#' @param R the number of bootstrap replications.
#' @param start `logical`; should `coef(x)` be passed as `start` to the `update(x, weights = ...)` call? In case the model `x` is computed by some numeric iteration, this may speed up the bootstrapping.
#' @param ... ignored.
#' @param fix `logical`; if `TRUE`, the covariance matrix is fixed to be positive semi-definite in case it is not.
#' @param use `character`; specification passed to [stats::cov()] for handling missing coefficients/parameters.
#'
#' @inherit sandwich::vcovBS return
#'
#' @details `vcovFWB()` functions like other `vcov()`-like functions, such as those in the `sandwich` package, in particular, \pkgfun{sandwich}{vcovBS}, which implements the traditional bootstrap (and a few other bootstrap varieties for linear models). Sets of weights are generated as described in the documentation for [fwb()], and the supplied model is re-fit using those weights. When the fitted model already has weights, these are multiplied by the bootstrap weights.
#'
#' For `lm` objects, the model is re-fit using [.lm.fit()] for speed, and, similarly, `glm` objects are re-fit using [glm.fit()] (or whichever fitting method was originally used). For other objects, [update()] is used to populate the weights and re-fit the model (this assumes the fitting function accepts non-integer case weights through a `weights` argument). If a model accepts weights in some other way, [fwb()] should be used instead; `vcovFWB()` is inherently limited in its ability to handle all possible models. It is important that the original model was not fit using frequency weights (i.e., weights that allow one row of data to represent multiple full, identical, individual units).
#'
#' See \pkgfun{sandwich}{vcovBS} and \pkgfun{sandwich}{vcovCL} for more information on clustering covariance matrices, and see [fwb()] for more information on how clusters work with the fractional weighted bootstrap. When clusters are specified, each cluster is given a bootstrap weight, and all members of the cluster are given that weight; estimation then proceeds as normal. By default, when `cluster` is unspecified, each unit is considered its own cluster.
#'
#' @seealso [fwb()] for performing the fractional weighted bootstrap on an arbitrary quantity; [fwb.ci()] for cmputing nonparametric confidence intervals for `fwb` objects; [summary.fwb()] for producing standard errors and confidence intervals for `fwb` objects; \pkgfun{sandwich}{vcovBS} for computing covariance matrices using the traditional bootstrap
#'
#' @export
#'
#' @examplesIf requireNamespace("lmtest", quietly = TRUE)
#' data("infert")
#' fit <- glm(case ~ spontaneous + induced, data = infert,
#'              family = "binomial")
#' lmtest::coeftest(fit, vcov. = vcovFWB, R = 200)
#'
#' @examplesIf requireNamespace("sandwich", quietly = TRUE)
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
vcovFWB <- function(x, cluster = NULL, R = 1000, start = FALSE, ..., fix = FALSE, use = "pairwise.complete.obs",
                    verbose = FALSE, cl = NULL) {

  #Check arguments
  chk::chk_count(R)
  chk::chk_flag(start)
  chk::chk_flag(fix)
  chk::chk_string(use)
  chk::chk_flag(verbose)

  ## set up return value with correct dimension and names
  cf <- coef(x)
  k <- length(cf)
  n <- nobs0(x)
  rval <- matrix(0, nrow = k, ncol = k, dimnames = list(names(cf), names(cf)))
  cf <- matrix(NA_real_, ncol = k, nrow = R, dimnames = list(NULL, names(cf)))

  ## cluster can either be supplied explicitly or
  ## be an attribute of the model...
  if (is.null(cluster)) cluster <- attr(x, "cluster")

  ## resort to cross-section if no clusters are supplied
  if (is.null(cluster)) cluster <- seq_len(n)

  ## collect 'cluster' variables in a data frame
  if (inherits(cluster, "formula")) {
    cluster_tmp <- suppressWarnings(expand.model.frame(x, cluster, na.expand = FALSE))
    cluster <- model.frame(cluster, cluster_tmp, na.action = na.pass)
  } else {
    cluster <- as.data.frame(cluster)
  }

  ## handle omitted or excluded observations
  if ((n != nrow(cluster)) && !is.null(x$na.action) && (class(x$na.action) %in% c("exclude", "omit"))) {
    cluster <- cluster[-x$na.action, , drop = FALSE]
  }

  if (nrow(cluster) != n) chk::err("number of observations in 'cluster' and 'nobs()' do not match")

  ## catch NAs in cluster -> need to be addressed in the model object by the user
  if (anyNA(cluster)) chk::err("cannot handle NAs in 'cluster': either refit the model without the NA observations in 'cluster' or impute the NAs")

  ## for multi-way clustering: set up interaction patterns
  p <- ncol(cluster)
  if (p > 1L) {
    clu <- lapply(seq_len(p), function(i) utils::combn(seq_len(p), i, simplify = FALSE))
    clu <- unlist(clu, recursive = FALSE)
    sign <- sapply(clu, function(i) (-1L)^(length(i) + 1L))
    paste_ <- function(...) paste(..., sep = "_")
    for (i in (p + 1L):length(clu)) {
      cluster <- cbind(cluster, Reduce(paste_, unclass(cluster[, clu[[i]] ]))) ## faster than: interaction()
    }
  } else {
    clu <- list(1)
    sign <- 1
  }

  opb <- pbapply::pboptions(type = if (verbose) "timer" else "none")
  on.exit(pbapply::pboptions(opb))

  ## apply infrastructure for refitting models
  applyfun <- function(X, FUN, ...) pbapply::pblapply(X, FUN, ..., cl = cl)

  ## use starting values?
  start <- if (start && inherits(x, "glm")) coef(x) else NULL

  ## bootstrap for each cluster dimension
  for (i in seq_along(clu)) {

    cli <- factor(cluster[[i]])

    ## bootstrap fitting function
    bootfit <- make.bootfit(x, cli, start)

    ## actually refit
    cf <- do.call("rbind", applyfun(seq_len(R), bootfit, ...))

    ## aggregate across cluster variables
    rval <- rval + sign[i] * cov(cf, use = use)
  }

  if (all_the_same(c(0, rval))) {
    chk::wrn("all variances and covariances are 0, indicating that the model failed to incorporate the bootstrapped weights")
    fix <- FALSE
  }

  ## check (and fix) if sandwich is not positive semi-definite
  if (fix && any((eig <- eigen(rval, symmetric = TRUE))$values < 0)) {
    eig$values <- pmax(eig$values, 0)
    rval[] <- crossprod(sqrt(eig$values) * t(eig$vectors))
  }
  return(rval)
}

nobs0 <- function (x, ...) {
  nobs1 <- stats::nobs
  nobs2 <- function(x, ...) NROW(residuals(x, ...))
  rval <- try(nobs1(x, ...), silent = TRUE)
  if (inherits(rval, "try-error") || is.null(rval))
    rval <- nobs2(x, ...)
  return(rval)
}

make.bootfit <- function(fit, cli, start, ...) {
  if (!is.factor(cli)) cli <- factor(cli)
  nc <- nlevels(cli)
  cluster_numeric <- as.integer(cli)

  if (identical(class(fit), "lm")) {
    mf <- model.frame(fit)
    x <- model.matrix(fit)
    y <- model.response(mf)
    if (!is.null(fit$offset)) y <- y - fit$offset

    bootfit <- function(j, ...) {
      #Generate cluster weights, assign to units
      cluster_w <- rexp(nc)
      cluster_w <- nc * cluster_w / sum(cluster_w)
      w <- cluster_w[cluster_numeric]

      if (!is.null(w0 <- weights(fit))) w <- w * w0

      ws <- sqrt(w)
      .lm.fit(x * ws, y = y * ws)$coefficients
    }
  }
  else if (inherits(fit, "glm")) {
    x <- model.matrix(fit)
    y <- fit[["y"]]

    if (is.null(fit[["method"]])) fit.fun <- stats::glm.fit
    else if (is.function(fit[["method"]])) fit.fun <- fit[["method"]]
    else if (is.character(fit[["method"]])) {
      fit.fun <- get0(fit[["method"]], envir = environment(fit[["terms"]]),
                      mode = "function")
      if (is.null(fit.fun)) {
        chk::err(sprintf("the `method` used to fit the original model (\"%s\") is unavailable", fit[["method"]]))
      }
    }
    else {
      chk::err("unrecognized fitting method; the model cannot be re-fit")
    }

    bootfit <- function(j, ...) {
      #Generate cluster weights, assign to units
      cluster_w <- rexp(nc)
      cluster_w <- nc * cluster_w / sum(cluster_w)
      w <- cluster_w[cluster_numeric]

      if (!is.null(w0 <- weights(fit))) w <- w * w0

      safe.glm.fit(fit.fun, x = x, y = y, weights = w, start = start,
                   offset = fit$offset, family = fit$family,
                   control = fit$control,
                   intercept = attr(fit$terms, "intercept") > 0)$coefficients
    }
  }
  else {
    bootfit <- function(j, ...) {
      #Generate cluster weights, assign to units
      cluster_w <- rexp(nc)
      cluster_w <- nc * cluster_w / sum(cluster_w)
      w <- cluster_w[cluster_numeric]

      if (!is.null(w0 <- weights(fit))) w <- w * w0

      up <- if (is.null(start)) {
        update(fit, weights = w, evaluate = TRUE)
      } else {
        update(fit, weights = w, start = start, evaluate = TRUE)
      }

      coef(up)
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
        !startsWith(conditionMessage(w), "non-integer x =")) warning(w)
    invokeRestart("muffleWarning")
  })
}
