expect_not_equal <- function(object, expected, ...,
                             tolerance = if (edition_get() >= 3) testthat_tolerance(),
                             info = NULL, label = NULL, expected.label = NULL) {

  if (!capabilities("long.double")) {
    return(NULL)
  }

  act <- quasi_label(rlang::enquo(object), label, arg = "object")
  exp <- quasi_label(rlang::enquo(expected), expected.label, arg = "expected")

  if (edition_get() >= 3) {
    expect_waldo_not_equal("equal", act, exp, info, ..., tolerance = tolerance)
  }
  else {
    if (!is.null(tolerance)) {
      comp <- compare(act$val, exp$val, ..., tolerance = tolerance)
    }
    else {
      comp <- compare(act$val, exp$val, ...)
    }
    expect(!comp$equal, sprintf("%s equal to %s.\n%s", act$lab, exp$lab, comp$message), info = info)
    invisible(act$val)
  }
}

expect_waldo_not_equal <- function(type, act, exp, info, ...) {
  comp <- waldo::compare(act$val, exp$val, ..., x_arg = "actual",
                         y_arg = "expected")

  expect(length(comp) > 0, sprintf("%s (%s) is %s to %s (%s).\n\n%s",
                                   act$lab, "`actual`", type, exp$lab,
                                   "`expected`",
                                   paste0(comp, collapse = "\n\n")),
         info = info, trace_env = rlang::caller_env())
  invisible(act$val)
}

#Use regex to make strings invariant to white spaces; also escape perens for perl
.w <- function(x) {
  x <- gsub("(", "\\(", x, fixed = TRUE)
  x <- gsub(")", "\\)", x, fixed = TRUE)
  gsub(" ", "(\\s+)", x, fixed = TRUE)
}

#The p-value column is named `Pr(>|t|)` when `ci.type = "cheap"` (the p-value comes
#from a t distribution) and `Pr(>|z|)` otherwise, so it has to be located rather than
#named. It is always the last column of `summary()`'s output.
p_value_column <- function(s) {
  s[, ncol(s)]
}

check_p_value_okay <- function(boot, ci.type, level, index, simultaneous = FALSE, eps) {
  R <- boot[["R"]]
  if (simultaneous) {

    eps_ <- max(eps, 1e-3)

    set.seed(7)
    suppressWarnings({
      s0 <- summary(boot, conf = level, ci.type = ci.type,
                   simultaneous = TRUE)
    })

    expect_s3_class(s0, "summary.fwb")
    expect_true(attr(s0, "simultaneous", TRUE))

    ## lower bound
    suppressWarnings({
      p <- p_value_column(summary(boot, conf = 0, p.value = TRUE,
                                  ci.type = ci.type, simultaneous = TRUE,
                                  null = s0[index, 3]))[index]
    })

    suppressWarnings({
      s <- summary(boot, conf = 1 - p, ci.type = ci.type,
                   simultaneous = TRUE)
    })

    expect_equal(s[index, 3:4], s0[index, 3:4],
                 tolerance = eps_,
                 ignore_attr = TRUE)

    ## upper bound
    suppressWarnings({
      p <- p_value_column(summary(boot, conf = 0, p.value = TRUE,
                                  ci.type = ci.type, simultaneous = TRUE,
                                  null = s0[index, 4]))[index]
    })

    suppressWarnings({
      s <- summary(boot, conf = 1 - p, ci.type = ci.type,
                   simultaneous = TRUE)
    })

    expect_equal(s[index, 3:4], s0[index, 3:4],
                 tolerance = eps_,
                 ignore_attr = TRUE)
  }
  else {
    suppressWarnings({
      s <- summary(boot, conf = level, ci.type = ci.type,
                   index = index)
    })

    expect_s3_class(s, "summary.fwb")
    expect_false(attr(s, "simultaneous", TRUE))

    ## lower bound
    suppressWarnings({
      p <- p_value_column(summary(boot, conf = 0, p.value = TRUE,
                                  ci.type = ci.type, index = index,
                                  null = s[1, 3]))[1]
    })

    expect_equal(p, 1 - level, tolerance = eps)

    ## upper bound
    suppressWarnings({
      p <- p_value_column(summary(boot, conf = 0, p.value = TRUE,
                                  ci.type = ci.type, index = index,
                                  null = s[1, 4]))[1]
    })

    expect_equal(p, 1 - level, tolerance = eps)
  }
}

# ---------------------------------------------------------------------------
# Condition expectations
#
# Messages raised by `arg::err()`/`arg::wrn()` are formatted by cli, which
# capitalizes the first letter, appends a period, converts inline markup (e.g.,
# `{.arg x}` to `` `x` ``), and hard-wraps the result to the console width. The
# wrapping means a literal `fixed = TRUE` match against a long message fails, and
# building a regex from the message means escaping every metacharacter cli may
# have introduced. `expect_err()`/`expect_wrn()` sidestep both: they collapse all
# whitespace in the *observed* message and then match a literal substring
# case-insensitively. Never copy a source string from `R/` into these -- use the
# rendered text.
# ---------------------------------------------------------------------------

#Collapse runs of whitespace so matching is invariant to how cli wrapped the message.
collapse_ws <- function(x) {
  gsub("\\s+", " ", trimws(paste(x, collapse = " ")))
}

#`cnds` is a list, because an expression can raise several warnings and the one under
#test is not always the first. The assertion is that *some* condition raised matched;
#a failure reports all of them, since which ones fired is the useful diagnostic.
.expect_cnd_text <- function(cnds, text, what) {
  #`as.list()` would decompose a lone condition object into its components, so callers
  #with a single condition are expected to wrap it themselves.
  cnds <- Filter(Negate(is.null), cnds)

  if (length(cnds) == 0L) {
    testthat::fail(sprintf("Expected %s, but none was signaled.", what))
    return(invisible(NULL))
  }

  if (is.null(text)) {
    testthat::succeed()
    return(invisible(cnds[[1L]]))
  }

  msgs <- vapply(cnds, function(cnd) collapse_ws(conditionMessage(cnd)), character(1L))

  testthat::expect_true(
    any(grepl(tolower(text), tolower(msgs), fixed = TRUE)),
    info = sprintf("no %s message contained the expected text.\n  expected: %s\n  actual:   %s",
                   what, encodeString(text, quote = "\""),
                   paste(encodeString(msgs, quote = "\""), collapse = "\n            "))
  )

  invisible(cnds[[1L]])
}

#Expect an error whose message contains `text` (literal, whitespace-insensitive).
expect_err <- function(object, text = NULL) {
  cnd <- NULL

  utils::capture.output(
    withCallingHandlers(cnd <- tryCatch({
      force(object)
      NULL
    }, error = function(e) e),
    message = function(m) invokeRestart("muffleMessage"),
    warning = function(w) invokeRestart("muffleWarning")))

  .expect_cnd_text(list(cnd), text, "error")
}

#Expect a warning whose message contains `text`. The expression still runs to
#completion, so assignments inside it take effect in the calling environment.
expect_wrn <- function(object, text = NULL) {
  cnds <- list()

  utils::capture.output(
    withCallingHandlers(force(object),
                        message = function(m) invokeRestart("muffleMessage"),
                        warning = function(w) {
                          cnds[[length(cnds) + 1L]] <<- w
                          invokeRestart("muffleWarning")
                        }))

  .expect_cnd_text(cnds, text, "warning")
}

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

#The tolerance used throughout. Bootstrap replicability comparisons are exact in
#principle -- the same weights fed to the same estimator -- so the only slack
#needed is for platforms without extended precision.
fwb_eps <- function() {
  if (capabilities("long.double")) 1e-8 else 1e-3
}

#A small deterministic dataset. Built without the RNG on purpose: a fixture that
#calls `set.seed()` would silently reset the stream that the replicability tests
#are about to pin down.
#
#`g` (cluster) is nested within `s` (stratum) because `g` cycles with a period
#that is a multiple of `s`'s -- `fwb()` requires that combination and rejects it
#otherwise.
fwb_test_df <- function(n = 60L) {
  x1 <- qnorm(ppoints(n))
  x2 <- as.numeric(seq_len(n) %% 3L == 0L)

  data.frame(x1 = x1,
             x2 = x2,
             y = 1 + x1 - 0.5 * x2 + sinpi(seq_len(n) / 7),
             yb = as.integer(x1 + sinpi(seq_len(n) / 5) > 0),
             s = factor(rep(c("a", "b", "c"), length.out = n)),
             g = factor(rep(seq_len(12L), length.out = n)))
}

#Statistics used by the tests. `w_identity_stat()` returns the bootstrap weights
#themselves, which turns "were the right weights used?" into an exact numeric
#comparison -- that is what makes the `fwb.array()` round-trip tests possible.
lm_stat <- function(data, w) {
  coef(lm(y ~ x1 + x2, data = data, weights = w))
}

glm_stat <- function(data, w) {
  coef(glm(yb ~ x1 + x2, data = data, family = quasibinomial(), weights = w))
}

w_identity_stat <- function(data, w) {
  as.numeric(w)
}

#A statistic with a random component, used to test the `random_statistic`
#detection and the code paths that depend on it. The perturbation is tiny so the
#estimates stay usable, but `identical()` still sees two different values.
random_lm_stat <- function(data, w) {
  coef(lm(y ~ x1 + x2, data = data, weights = w)) + rnorm(3L, sd = 1e-7)
}

#Reparent the statistics to the global environment.
#
#`statistic` is serialized to parallel workers along with its enclosing environment.
#testthat sources this file into an environment whose parent is the *attached*
#`package:fwb`, and base R warns when it serializes a reference to an attached package
#environment ("'package:fwb' may not be available when loading") -- once per worker per
#call, which buried the real output of the backend tests in warnings.
#
#Nothing above needs anything but base and stats, and a user writing one of these at the
#top level of a script would get `globalenv()` anyway, so this also makes the tests
#reflect ordinary use. Note that it is the attached package environment that triggers
#this, not a namespace: a statistic defined inside a package does not warn.
for (.stat in c("lm_stat", "glm_stat", "w_identity_stat", "random_lm_stat")) {
  assign(.stat, `environment<-`(get(.stat), globalenv()))
}

rm(.stat)

# ---------------------------------------------------------------------------
# Parallel-backend scaffolding
# ---------------------------------------------------------------------------

#The largest number of parallel workers the tests may ask for.
#
#`R CMD check --as-cran` sets `_R_CHECK_LIMIT_CORES_`, and `parallel:::.check_ncores()`
#then *errors* on a request for more than two processes -- CRAN policy is a hard maximum
#of two. `devtools::check()` also sets `NOT_CRAN = "true"`, so `skip_on_cran()` does not
#skip and the backend tests run under that limit: asking for three there is a failure
#rather than a skip. That is the whole reason `devtools::test()` can be clean while
#`devtools::check()` is not, so it cannot be diagnosed from `devtools::test()` alone.
#
#The condition mirrors `.check_ncores()` itself: the variable counts as set unless it is
#literally `"false"`, and `"warn"` still means the limit applies.
#
#Capping rather than skipping keeps the "the results do not depend on the worker count"
#comparison running wherever a third worker is allowed, and degrades it to a two-worker
#comparison where one is not. Two distinct counts are all that comparison needs; the
#third is only ever a second data point. Callers build their labels from the value so a
#failure never reports a worker count that was not actually used.
max_test_workers <- function() {
  egc <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

  if (nzchar(egc) && !identical(egc, "false")) {
    return(2L)
  }

  3L
}

#Every backend helper below raises the `parallelly` soft limit on workers.
#Without it, CI machines that report a single core turn each `plan()`/
#`makeCluster()` call into a warning, and the `expect_no_condition()`
#expectations in these tests would fail on the scaffolding rather than on `fwb()`.
local_worker_limit <- function(.env = parent.frame()) {
  rlang::local_options(parallelly.maxWorkers.localhost = Inf,
                       parallelly.availableCores.min = 2L,
                       .frame = .env)
}

#Register a cleanup call as an `on.exit()` expression in the calling frame.
#`rlang::defer()` is not exported and `withr` is not a declared dependency, so
#this keeps the scaffolding dependency-free.
#
#`template` is a quoted expression in which `.` stands for the object to clean up;
#`substituteDirect()` splices that object in as a *value*, so nothing depends on
#what the caller named it or on the binding still existing at exit time.
defer_call <- function(template, x, .env = parent.frame()) {
  expr <- substituteDirect(template, list(. = x))

  do.call(base::on.exit,
          list(expr, add = TRUE, after = FALSE),
          envir = .env)
}

#Set a `future` plan for the duration of the calling test and restore the
#previous one afterwards. Environments that cannot start workers (no sockets, no
#forking) skip rather than fail -- these tests are about `fwb()`'s use of the
#backend, not about whether the backend can start.
local_plan <- function(strategy = "multisession", workers = 2L, .env = parent.frame()) {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")

  local_worker_limit(.env = .env)

  old <- rlang::try_fetch(
    suppressWarnings({
      if (identical(strategy, "sequential")) future::plan(future::sequential)
      else future::plan(strategy, workers = workers)
    }),
    error = function(e) {
      skip(sprintf("could not set up a %s future plan: %s", strategy,
                   conditionMessage(e)))
    }
  )

  defer_call(quote(suppressWarnings(future::plan(.))), old, .env = .env)

  invisible(old)
}

#As above, for a PSOCK cluster object.
local_cluster <- function(workers = 2L, .env = parent.frame()) {
  skip_if_not_installed("parallel")

  local_worker_limit(.env = .env)

  cl <- rlang::try_fetch(
    suppressWarnings(parallel::makePSOCKcluster(workers)),
    error = function(e) {
      skip(sprintf("could not start a PSOCK cluster: %s", conditionMessage(e)))
    }
  )

  #`try()` so the cleanup is idempotent: a test that needs to show behavior after the
  #workers are gone stops the cluster itself, and stopping an already-stopped cluster
  #errors with "invalid connection".
  defer_call(quote(try(parallel::stopCluster(.), silent = TRUE)), cl, .env = .env)

  cl
}

#`cl = <integer>` routes to `parallel::mclapply()`, which forks.
skip_if_no_forking <- function() {
  if (.Platform$OS.type == "windows") {
    skip("`cl = <integer>` falls back to sequential evaluation on Windows")
  }

  invisible(NULL)
}

# ---------------------------------------------------------------------------
# progressr
# ---------------------------------------------------------------------------

#Count the progress updates signaled while `expr` runs.
#
#`enable = TRUE` is not optional: `progressr` reports progress only when
#`interactive()` by default, so without it every count comes back 0 during a test
#run and the assertions below would pass no matter what `progressify()` did.
count_progress_updates <- function(expr) {
  skip_if_not_installed("progressr")

  n <- 0L

  progressr::with_progress(
    withCallingHandlers(force(expr),
                        condition = function(cnd) {
                          if (inherits(cnd, "progression") &&
                              identical(cnd[["type"]], "update")) {
                            n <<- n + 1L
                          }
                        }),
    handlers = progressr::handler_void(),
    enable = TRUE)

  n
}

# ---------------------------------------------------------------------------
# Comparisons
# ---------------------------------------------------------------------------

#Compare the bootstrap replicates of two `<fwb>` objects. Names are dropped
#because `t`'s dimnames come from `t0`, which is not what any of these tests are
#about.
expect_same_t <- function(x, y, tolerance = fwb_eps(), ...) {
  expect_equal(unname(x[["t"]]), unname(y[["t"]]),
               tolerance = tolerance, ...)
}

expect_different_t <- function(x, y, tolerance = fwb_eps(), ...) {
  expect_not_equal(unname(x[["t"]]), unname(y[["t"]]),
                   tolerance = tolerance, ...)
}

#Send plot output to a null device for the duration of the calling test. `setup.R`
#opens one for the whole run, but a test that inspects `par()` needs a device it
#controls, and one that errors mid-plot should not leave the shared device in a
#modified state.
local_null_device <- function(.env = parent.frame()) {
  grDevices::pdf(NULL)

  do.call(base::on.exit,
          list(quote(grDevices::dev.off()), add = TRUE, after = FALSE),
          envir = .env)

  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Simultaneous inference
# ---------------------------------------------------------------------------

#Empirical simultaneous coverage of the pointwise `ci.type` intervals at level `l`: the
#proportion of bootstrap draws whose entire parameter vector falls inside them.
#
#This lives here rather than in the package because it is the *independent* route to a
#quantity the package computes a much faster way. `simultaneous_ci_level()` never builds an
#interval -- it reads the answer off the order statistics -- so checking it against a
#definition written in terms of `compute_ci()` output is what tests the algebra rather than
#restating it.
simultaneous_coverage <- function(object, index, level, ci.type = "perc") {
  estimates <- t(object[["t"]][, index, drop = FALSE])
  k <- length(index)

  suppressWarnings({
    ci.out <- fwb:::compute_ci(ci.type, t = object[["t"]], t0 = object[["t0"]],
                               conf = level, index = index, boot.out = object)
  })

  nc <- ncol(ci.out)

  interval <- ci.out[, c(nc - 1L, nc), drop = FALSE]

  mean(colSums(estimates >= interval[, 1L]) == k &
         colSums(estimates <= interval[, 2L]) == k)
}
