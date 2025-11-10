#' Fractional Weighted Bootstrap Covariance Matrix Estimation
#'
#' `vcovFWB()` estimates the covariance matrix of model coefficient estimates using the fractional weighted bootstrap. It serves as a drop-in for `[stats::vcov()]` or `sandwich::vcovBS()`. Clustered covariances are can be requested.
#'
#' @inheritParams sandwich::vcovBS
#' @inheritParams fwb
#' @param x a fitted model object, such as the output of a call to `lm()` or `glm()`. The model object must result from a function that can be updated using [update()] and has a `weights` argument to input non-integer case weights.
#' @param R the number of bootstrap replications. Default is 1000 (more is better but slower).
#' @param start `logical`; should `.coef(x)` be passed as `start` to the `update(x, weights = ...)` call? In case the model `x` is computed by some numeric iteration, this may speed up the bootstrapping. Default is `FALSE`.
#' @param wtype string; the type of weights to use. Allowable options include `"exp"` (the default), `"pois"`, `"multinom"`, `"mammen"`, `"beta"`, and `"power"`. See [fwb()] for details. See [set_fwb_wtype()] to set a global default.
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
#' @seealso
#' * [fwb()] for performing the fractional weighted bootstrap on an arbitrary quantity
#' * [fwb.ci()] for computing nonparametric confidence intervals for `fwb` objects
#' * [summary.fwb()] for producing standard errors and confidence intervals for `fwb` objects
#' * \pkgfun{sandwich}{vcovBS} for computing covariance matrices using the traditional bootstrap (the fractional weighted bootstrap is also available but with limited options).
#'
#' @examplesIf rlang::is_installed("lmtest")
#' set.seed(123, "L'Ecuyer-CMRG")
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
  chk::chk_count(R)
  chk::chk_flag(start)
  chk::chk_flag(fix)
  chk::chk_string(use)
  chk::chk_flag(verbose)
  chk::chk_function(.coef)

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

  ## set up return value with correct dimension and names
  cf <- .coef(x)

  if (!is.numeric(cf) || length(dim(cf)) > 1L) {
    if (identical(.coef, eval(formals()[[".coef"]]))) {
      .err("the coefficients extracted using `coef()` from the supplied model are not in the form of a numeric vector; see the `.coef` argument at `help(\"vcovFWB\")`")
    }
    else {
      .err("the function supplied to `.coef` must return a numeric vector")
    }
  }

  k <- length(cf)
  n <- nobs0(x)
  rval <- matrix(0, nrow = k, ncol = k, dimnames = list(names(cf), names(cf)))
  cf <- matrix(NA_real_, ncol = k, nrow = R, dimnames = list(NULL, names(cf)))

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
  if (n != nrow(cluster) && is_not_null(x$na.action) && inherits(x$na.action, c("exclude", "omit"))) {
    cluster <- cluster[-x$na.action, , drop = FALSE]
  }

  if (nrow(cluster) != n) {
    .err("number of observations in `cluster` and `nobs()` do not match")
  }

  ## catch NAs in cluster -> need to be addressed in the model object by the user
  if (anyNA(cluster)) {
    .err("cannot handle `NA`s in `cluster`: either refit the model without the `NA` observations in `cluster` or impute the `NA`s")
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

  opb <- pbapply::pboptions(type = if (verbose) "timer" else "none")
  on.exit(pbapply::pboptions(opb))

  .env <- parent.frame()

  ## apply infrastructure for refitting models
  applyfun <- {
    if (identical(cl, "future"))
      function(X, FUN, ...) pbapply::pblapply(X, FUN, ..., cl = cl, future.seed = TRUE)
    else
      function(X, FUN, ...) pbapply::pblapply(X, FUN, ..., cl = cl)
  }

  ## use starting values?
  start <- if (start && inherits(x, "glm")) .coef(x) else NULL

  ## bootstrap for each cluster dimension
  for (i in seq_along(clu)) {

    cli <- factor(cluster[[i]])

    ## bootstrap fitting function
    bootfit <- make.bootfit(x, cli, start, drop0 = drop0,
                            gen_weights, .coef, .env)

    ## actually refit
    cf <- do.call("rbind", applyfun(seq_len(R), bootfit, ...))

    ## aggregate across cluster variables
    rval <- rval + sgn[i] * stats::cov(cf, use = use)
  }

  if (all_the_same(c(0, rval))) {
    .wrn("all variances and covariances are 0, indicating that the model failed to incorporate the bootstrapped weights")
    fix <- FALSE
  }

  ## check (and fix) if sandwich is not positive semi-definite
  if (fix) {
    eig <- eigen(rval, symmetric = TRUE)$values

    if (any(eig < 0)) {
      eig$values <- pmax(eig$values, 0)
      rval[] <- crossprod(sqrt(eig$values) * t(eig$vectors))
    }
  }

  rval
}

nobs0 <- function(x, ...) {
  if (inherits(x, "coxph")) {
    rval <- x[["n"]]
  }
  else {
    rval <- try(stats::nobs(x, ...), silent = TRUE)
  }

  if (is_null(rval) || inherits(rval, "try-error")) {
    rval <- NROW(residuals(x, ...))
  }

  rval
}

make.bootfit <- function(fit, cli, start, drop0, gen_weights, .coef, .env) {
  cli <- as.factor(cli)
  nc <- nlevels(cli)
  cluster_numeric <- as.integer(cli)
  special_coef <- !identical(.coef(fit), try(stats::coef(fit), silent = TRUE))
  bootfit <- NULL

  w0 <- weights(fit) %or% 1

  if (!special_coef) {
    # Test if model uses w_*() functions by seeing if weighted
    # model.matrix differs from original; if so, need to use
    # general update instead of lm- or glm-specific.

    mm0 <- model.matrix(fit)

    if (is_not_null(mm0)) {
      bootfit_test <- function(j, ...) {
        rlang::local_options(fwb_internal_w_env = rlang::current_env())

        #Generate cluster weights, assign to units
        cluster_w <- drop(gen_weights(nc, 1L))
        .wi <- cluster_w[cluster_numeric]

        .wi <- .set_class(.wi, "fwb_internal_w")

        w <- .wi * w0

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

        up <- do.call("update", args)

        utils::capture.output({
          up <- eval(up, envir = .env)
        })

        mm <- try(model.matrix(up), silent = TRUE)

        if (inherits(mm, "try-error")) {
          return(NULL)
        }

        mm
      }

      with_seed_preserved({
        mm_test <- bootfit_test(1L)
      })

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

    bootfit <- function(j, ...) {
      #Generate cluster weights, assign to units
      cluster_w <- drop(gen_weights(nc, 1L))
      .wi <- cluster_w[cluster_numeric]

      .wi <- .wi * w0

      ws <- sqrt(.wi)

      .lm.fit(x * ws, y = y * ws)$coefficients
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
        .err(sprintf("the `method` used to fit the original model (%s) is unavailable",
                     add_quotes(fit[["method"]])))
      }
    }
    else {
      .err("unrecognized fitting method; the model cannot be re-fit")
    }

    bootfit <- function(j, ...) {
      #Generate cluster weights, assign to units
      cluster_w <- drop(gen_weights(nc, 1L))
      .wi <- cluster_w[cluster_numeric]

      .wi <- .wi * w0

      if (!isFALSE(drop0)) {
        zero_w <- .wi == 0

        if (any(zero_w)) {
          if (isTRUE(drop0)) {
            .wi <- .wi[!zero_w]

            x <- x[!zero_w, , drop = FALSE]
            y <- y[!zero_w]
            offset <- offset[!zero_w]
          }
          else {
            is.na(.wi)[!zero_w] <- TRUE
          }
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
      cluster_w <- drop(gen_weights(nc, 1L))
      .wi <- cluster_w[cluster_numeric]

      .wi <- .set_class(.wi, "fwb_internal_w")

      w <- .wi * w0

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
      .wrn(w, tidy = FALSE, immediate = FALSE)
    }
    invokeRestart("muffleWarning")
  })
}
