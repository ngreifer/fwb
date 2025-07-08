#' @title Bearing Cage field failure data
#' @description The data consist of 1703 aircraft engines put into service over time. There were 6 failures and 1697 right-censored observations. These data were originally given in Abernethy et al. (1983) and were reanalyzed in Meeker and Escobar (1998, chap.8). The dataset used here specifically comes from Xu et al. (2020) and is used in a Weibull analysis of failure times.
#' @usage data("bearingcage")
#' @format A data frame with 1703 rows and 2 variables:
#' \describe{
#'   \item{`hours`}{`integer`; the number of hours until failure or censoring}
#'   \item{`failure`}{`logical`; whether a failure occurred}
#'}
#' @references
#' Abernethy, R. B., Breneman, J. E., Medlin, C. H., and Reinman, G. L. (1983), "Weibull Analysis Handbook," Technical Report, Air Force Wright Aeronautical Laboratories. \doi{10.21236/ADA143100}
#'
#' Meeker, W. Q., and Escobar, L. A. (1998), *Statistical Methods for Reliability Data*, New York: Wiley.
#'
#' Xu, L., Gotwalt, C., Hong, Y., King, C. B., & Meeker, W. Q. (2020). Applications of the Fractional-Random-Weight Bootstrap. *The American Statistician*, 74(4), 345–358. \doi{10.1080/00031305.2020.1731599}
#'
"bearingcage"
