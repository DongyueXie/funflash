# Provides functions to solve the Empirical Bayes Normal Means problem.
# Each function must take arguments x, s, ebnm_param, and output, and
# must return a list with elements postmean, postmean2, fitted_g,
# and penloglik.
#
# If sampling from the posterior is desired, then the function must also
# be able to return a sampling function when output = "post_sampler".
# This sampling function should take a single argument nsamp and return
# a matrix with nsamp rows and length(x) columns. (In other words, the
# (i,j)-entry of the matrix should correspond to the ith sample from the
# posterior for the jth element in the vector of observations.)


# @title EBNM using ash
#
# @description A wrapper to the ash function for flash.
#
# @param x A vector of observations.
#
# @param s A vector of standard errors.
#
# @param ash_param A list of parameters to be passed into ash.
#
# @param output If output = "post_sampler", then the return value is a
#   function that samples from the posterior.
#
#' @importFrom ashr ash
#'
ebnm_ash = function(x, s, ash_param, output = NULL) {
  if (identical(output, "post_sampler")) {
    ash_param$output = "post_sampler"
  } else {
    ash_param$output = "flash_data"
  }

  a = do.call(ash,
              c(list(betahat = as.vector(x), sebetahat = as.vector(s)),
                ash_param))

  if (identical(output, "post_sampler")) {
    out = a$post_sampler
  } else if (is.null(a$flash_data$postmean)) {
    stop(paste("ashr is not outputting flashr data in the right format.",
               "Maybe ashr needs updating to latest version?"))
  } else {
    out = a$flash_data
  }
  return(out)
}


# @title EBNM using point-exponential prior, from ebnm package
#
# @description A wrapper to the function
#   \code{\link[ebnm]{ebnm_point_exponential}}.
#
# @inheritParams ebnm_ash
#
# @param ebnm_param A list of parameters to be passed to
#   \code{ebnm_point_exponential}.
#
# Do not include an @importFrom field here until ebnm is on CRAN.
#
ebnm_pe = function(x, s, ebnm_param, output = NULL) {
  if (identical(output, "post_sampler")) {
    ebnm_param$output = "posterior_sampler"
  } else {
    ebnm_param$output = c("posterior_mean", "posterior_second_moment",
                          "fitted_g", "log_likelihood")
  }

  if (!is.null(ebnm_param$g)) {
    ebnm_param$g_init = ebnm_param$g
    ebnm_param$g = NULL
  }
  if (!is.null(ebnm_param$fixg)) {
    ebnm_param$fix_g = ebnm_param$fixg
    ebnm_param$fixg = NULL
  }

  res = do.call(ebnm::ebnm_point_exponential,
                c(list(x = as.vector(x), s = as.vector(s)),
                  ebnm_param))

  if (identical(output, "post_sampler")) {
    out = res$posterior_sampler
  } else {
    out = list(postmean = res$posterior$mean,
               postmean2 = res$posterior$second_moment,
               fitted_g = res$fitted_g,
               penloglik = res$log_likelihood)
  }

  return(out)
}

# @title EBNM using point-normal prior, from ebnm package
#
# @description A wrapper to the function
#   \code{\link[ebnm]{ebnm_point_normal}}.
#
# @inheritParams ebnm_ash
#
# @param ebnm_param A list of parameters to be passed to
#   \code{ebnm_point_normal}.
#
# Do not include an @importFrom field here until ebnm is on CRAN.
#
ebnm_pn = function(x, s, ebnm_param, output = NULL) {
  if (identical(output, "post_sampler")) {
    ebnm_param$output = "posterior_sampler"
  } else {
    ebnm_param$output = c("posterior_mean", "posterior_second_moment",
                          "fitted_g", "log_likelihood")
  }

  if (!is.null(ebnm_param$g)) {
    ebnm_param$g_init = ebnm_param$g
    ebnm_param$g = NULL
  }
  if (!is.null(ebnm_param$fixg)) {
    ebnm_param$fix_g = ebnm_param$fixg
    ebnm_param$fixg = NULL
  }

  res = do.call(ebnm::ebnm_point_normal,
                c(list(x = as.vector(x), s = as.vector(s)),
                  ebnm_param))

  if (identical(output, "post_sampler")) {
    out = res$posterior_sampler
  } else {
    out = list(postmean = res$posterior$mean,
               postmean2 = res$posterior$second_moment,
               fitted_g = res$fitted_g,
               penloglik = res$log_likelihood)
  }

  return(out)
}


# @title EBNM using point-laplace prior, from ebnm package
#
# @description A wrapper to the function
# \code{\link[ebnm]{ebnm_point_laplace}}.
#
# @inheritParams ebnm_pn
#
# @param ebnm_param A list of parameters to be passed to the function
#   \code{ebnm_point_laplace}.
#
# @param output Sampling from the posterior has not yet been implemented
#   for point-laplace priors.
#
# Do not include an @importFrom field here until ebnm is on CRAN.
#
ebnm_pl = function(x, s, ebnm_param, output = NULL) {
  if (identical(output, "post_sampler")) {
    ebnm_param$output = "posterior_sampler"
  } else {
    ebnm_param$output = c("posterior_mean", "posterior_second_moment",
                          "fitted_g", "log_likelihood")
  }

  if (!is.null(ebnm_param$g)) {
    ebnm_param$g_init = ebnm_param$g
    ebnm_param$g = NULL
  }
  if (!is.null(ebnm_param$fixg)) {
    ebnm_param$fix_g = ebnm_param$fixg
    ebnm_param$fixg = NULL
  }

  res = do.call(ebnm::ebnm_point_laplace,
                c(list(x = as.vector(x), s = as.vector(s)),
                  ebnm_param))

  return(list(postmean = res$posterior$mean,
              postmean2 = res$posterior$second_moment,
              fitted_g = res$fitted_g,
              penloglik = res$log_likelihood))
}


ebnm_dwt_haaar = function(x, s, g_init=NULL, fix_g=F, output){
  ebpmf:::ebnm_dwt(x, s, g_init, fix_g,filter.number=1,family="DaubExPhase")
}

# @title EBNM using wavelet prior
ebnm_smooth = function (x, s, ebnm_param, output = NULL)
{
  if (identical(output, "post_sampler")) {
    ebnm_param$output = "posterior_sampler"
  }
  else {
    ebnm_param$output = c("posterior_mean", "posterior_second_moment",
                          "fitted_g", "log_likelihood")
  }
  if (!is.null(ebnm_param$g)) {
    ebnm_param$g_init = ebnm_param$g
    ebnm_param$g = NULL
  }
  if (!is.null(ebnm_param$fixg)) {
    ebnm_param$fix_g = ebnm_param$fixg
    ebnm_param$fixg = NULL
  }
  res = do.call(ebnm_dwt_haaar, c(list(x = as.vector(x),
                                              s = as.vector(s)), ebnm_param))
  if (identical(output, "post_sampler")) {
    out = res$posterior_sampler
  }
  else {
    out = list(postmean = res$posterior$mean, postmean2 = res$posterior$second_moment,
               fitted_g = res$fitted_g, penloglik = res$log_likelihood)
  }
  return(out)
}



# @title EBNM using wavelet prior
ebnm_smooth_ndwt = function (x, s, ebnm_param, output = NULL)
{
  if (identical(output, "post_sampler")) {
    ebnm_param$output = "posterior_sampler"
  }
  else {
    ebnm_param$output = c("posterior_mean", "posterior_second_moment",
                          "fitted_g", "log_likelihood")
  }
  if (!is.null(ebnm_param$g)) {
    ebnm_param$g_init = ebnm_param$g
    ebnm_param$g = NULL
  }
  if (!is.null(ebnm_param$fixg)) {
    ebnm_param$fix_g = ebnm_param$fixg
    ebnm_param$fixg = NULL
  }
  res = do.call(ebpmf::ebnm_ndwt, c(list(x = as.vector(x),
                                       s = as.vector(s)), ebnm_param))
  if (identical(output, "post_sampler")) {
    out = res$posterior_sampler
  }
  else {
    out = list(postmean = res$posterior$mean, postmean2 = res$posterior$second_moment,
               fitted_g = res$fitted_g, penloglik = res$log_likelihood)
  }
  return(out)
}
