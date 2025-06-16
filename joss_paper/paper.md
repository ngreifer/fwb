---
title: 'fwb: An R package for the fractional weighted bootstrap'
tags:
  - R
authors:
  - name: Noah Greifer
    orcid: 0000-0003-3067-7154
    corresponding: true
    affiliation: 1
affiliations:
 - name: Institute for Quantitative Social Science, Harvard University, United States
   index: 1
   ror: 03vek6s52
date: 16 June 2025
bibliography: paper.bib
---

# Summary

The nonparametric bootstrap is a simple method for performing inference on estimates of statistical parameters without making strong assumptions about the data or a statistical model and without needing to derive the formula for an estimate's variance [@efronIntroductionBootstrap1993]. The traditional nonparametric bootstrap involves resampling the observations from the sample with replacement, computing the estimates in the bootstrap sample, performing these steps many times, and summarizing the resulting distribution of bootstrap estimates to compute their confidence intervals. This is equivalent to drawing integer weights from a multinomial distribution and computing weighted estimates in each bootstrap sample. However, integer weights are not the only weights that can be used, and fractional (non-integer) weights carry certain advantages, described below. `fwb` implements the fractional weighted bootstrap, also known as the Bayesian bootstrap, using a simple interface and with options for parallelization.

# Statement of need

The fractional weighted bootstrap was described originally by @rubinBayesianBootstrap1981 and, for a more applied audience, by @xuApplicationsFractionalRandomWeightBootstrap2020a. In each bootstrap sample, weights are drawn from a specified distribution, and the statistics of interest are computed in the weighted sample. In contrast to the traditional bootstrap, which drops any units that are not included in a given bootstrap sample in computing that sample's estimates, the fractional weighted bootstrap retains all units in every bootstrap sample. This yields several benefits, including that rare values are included in every bootstrap sample, and the resulting bootstrap distribution will be smoother (when the weights are drawn from a continuous distribution) [@xuApplicationsFractionalRandomWeightBootstrap2020a].

This is demonstrated in @xuApplicationsFractionalRandomWeightBootstrap2020a, who use the fractional weighted bootstrap to estimate the center and shape of a Weibull distribution characterizing failure times for a set of 1703 aircraft engines, among which only 6 failures were observed, data previously analyzed by @meekerStatisticalMethodsReliability1998. In cases with rare outcomes, rare covariate profiles, or highly influential observations, the traditional bootstrap can fail to yield valid inference. The fractional weighted bootstrap avoids these issues by retaining all units and allowing the weights to take on strictly positive, non-integer values.

The traditional bootstrap is implemented in the `boot` package [@canty2024; @davisonBootstrapMethodsTheir1997], which uses a simple interface to perform bootstrapping. The user supplies a dataset and a function that computes the desired statistics on each bootstrap sample from the dataset to `boot::boot()`. `boot::boot.ci()` is used to compute confidence intervals using a variety of methods. `fwb` provides a drop-in replacement for `boot::boot()` in `fwb::fwb()`, which similarly takes in the original dataset and a function that computes the desired statistics incorporating the bootstrap weights. `fwb()` supports parallelizaton through the `pbapply` package [@solymos2023], which provides interfaces to the `parallel` and `future` packages [@bengtsson2021] for parallel processing, as well as stratified bootstrapping (i.e., bootstrapping within strata), clustered bootstrapping (bootstrapping of clusters), and different distributions from which the weights are drawn. These distributions include the Dirichlet/Exponential weights originally described by @rubinBayesianBootstrap1981, "Mammen" weights that mimic the form of the distribution of weights proposed by @mammenBootstrapWildBootstrap1993 for the wild bootstrap, multinomial weights to exactly reproduce the traditional nonparametric bootstrap by `boot`, and Poisson frequency weights [@hanley2006].

Like `boot`, `fwb` provides tools for plotting the bootstrap distributions and computing confidence intervals, including the Normal approximation, percentile, bias-corrected percentile, and bias-corrected and accelerated (BCa) confidence intervals implemented in `fwb.ci()` [@davisonBootstrapMethodsTheir1997]. `fwb` also provides additional S3 methods for `fwb` objects, including a `summary()` method to facilitate extraction of estimates and confidence intervals and perform hypothesis tests, which are otherwise fairly challenging in `boot`.

A common use of the bootstrap is to compute an estimate of the asymptotic covariance matrix of parameter estimates in regression models. Bootstrapped confidence intervals avoid some of the often strict and unrealistic assumptions of maximum likelihood-based intervals, in particular, that the model is correctly specified (e.g., that the conditional variance is constant in a linear regression model) [@berkAssumptionLeanRegression2021]. The `sandwich` package [@zeileis2020] provides `sandwich::vcovBS()` to compute bootstrapped covariance matrices of the regression parameters, and `fwb` provides `vcovFWB()` for the same purpose but supporting the fractional weighted bootstrap[^1]. These functions can be used as a drop-in for `vcov()` for many regression functions (in particular, those that have an `update()` method). `fwb` also includes utility functions for computing weighted statistics and variable transformations, such as `w_mean()` for computing the weighted mean and `w_std()` for standardizing (scaling and centering), which automatically incorporate bootstrap weights when used inside `fwb()` or `vcovFWB()`.

[^1]: Note that `sandwich::vcovBS()` supports a basic form of the fractional weighted bootstrap but does not support different distributions of weights or other options provided by `fwb`.

The `bayesboot` package [@bayesboot] also implements the Bayesian bootstrap, but it does so in an explicitly Bayesian fashion, focusing on visualizing the bootstrap distribution as a Bayesian posterior distribution. In contrast, `fwb` is designed to be part of a frequentist workflow similar to traditional bootstrap implementations, e.g., in how it reports point estimates, confidence intervals, and hypothesis tests rather than using the bootstrap distribution to characterize the posterior of the estimate. `bayesboot` supports bootstrapping with Dirichlet weights and supports the use of statistics that cannot directly incorporate weights, but it does not support other weight distributions, clustering, or confidence interval estimation (other than highest density credible intervals that are part of standard Bayesian analysis workflows). `fwb` is designed to provide an alternative to `boot` for frequentist inference using the more robust fractional weighted bootstrap rather than implementing the Bayesian bootstrap as a specific Bayesian procedure.

`fwb` is available on [CRAN](https://cran.r-project.org/package=fwb).

# References
