#' Set weights type
#'
#' Set the default for the type of weights used in the weighted bootstrap computed by [fwb()] and [vcovFWB()].
#'
#' @param wtype string; the type of weights to use. Allowable options include `"exp"` (the default), `"pois"`, `"multinom"`, and `"mammen"`. Abbreviations allowed. See Details for what these mean.
#' @param fwb optional; an `fwb` objct, the output of a call to [fwb()]. If left empty, will extract the weights type from `options()`.
#'
#' @return `set_fwb_wtype()` returns a call to [options()] with the options set to those prior to `set_fwb_wtype()` being called. This makes it so that calling `options(op)`, where `op` is the output of `set_fwb_wtype()`, resets the `fwb_wtype` to its original value. `get_fwb_wtype()` returns a string containing the `fwb_wtype` value set globally (if no argument is supplied) or used in the supplied `fwb` object.
#'
#' @details `set_fwb_wtype(x)` is equivalent to calling `options(fwb_wtype = x)`. `get_fwb_wtype()` is equivalent to calling `getOption("fwb_wtype")` when no argument is supplied and to extracting the `wtype` component of an `fwb` object when supplied.
#'
#' @seealso [fwb] for a definition of each types of weights; [vcovFWB()]; [options()]; \pkgfun{boot}{boot} for the traditional bootstrap.
#'
#' @examplesIf requireNamespace("survival", quietly = TRUE)
#' # Performing a Weibull analysis of the Bearing Cage
#' # failure data as done in Xu et al. (2020)
#' set.seed(123, "L'Ecuyer-CMRG")
#' data("bearingcage")
#'
#' #Set fwb type to "mammen"
#' op <- set_fwb_wtype("mammen")
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
#' #Get the fwb type used in the bootstrap
#' get_fwb_wtype(boot_est)
#' get_fwb_wtype()
#'
#' #Restore original options
#' options(op)
#'
#' get_fwb_wtype()
#'

#' @export
set_fwb_wtype <- function(wtype = getOption("fwb_wtype", "exp")) {
  chk::chk_string(wtype)

  wtype <- tolower(wtype)
  wtype <- match_arg(wtype, c("exp", "multinom", "poisson", "mammen"))

  op <- options("fwb_wtype" = wtype)

  invisible(op)
}

#' @export
#' @rdname set_fwb_wtype
get_fwb_wtype <- function(fwb) {
  if (missing(fwb)) {
    return(getOption("fwb_wtype", "exp"))
  }

  if (!inherits(fwb, "fwb")) {
    .err("the argument to `get_fwb_wtype()` must either be left empty or the output of a call to `fwb()`")
  }

  fwb[["wtype"]]
}
