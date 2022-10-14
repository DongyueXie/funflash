#'@title Fit smoothed flash
#'@description This model fits wavelet prior for each factors.
#'@param Y data matrix, n by p, p is a power of 2.
#'@param S An n by p matrix of the standard errors of the observations in Y. (Can be a scalar if all standard errors are equal.) If S = NULL, then the standard errors will be estimated during fitting.
#'@param Kmax maximum number of topics
#'@param tol tolerance
#'@param init_fn initialization of
#'@param type wavelet or station. Use wavelet for exact ELBO calc.
#'@param var.type Use 'by_row' for exact ELBO calc
#'@return same as flashr
#'@import flashr
#'@import smashr
#'@export
smooth_flash = function(Y,S=NULL,Kmax=10,tol=0.01,
                        init_fn = "udv_si",
                        ebnm_fn = "ebnm_pn",
                        var_type = 'by_row',
                        ebnm_param = NULL,
                        f_init = NULL,
                        filter.number = 1,
                        family = 'DaubExPhase',
                        nullcheck=TRUE,
                        maxiter=100,
                        verbose = TRUE,
                        type='wavelet',
                        seed=12345){

  set.seed(seed)
  f = flashr:::handle_f(f_init, init_null_f = TRUE)
  data = flashr:::flash_set_data(Y,S)
  data = flashr:::handle_data(data,f)
  if(!is.null(S)){
    f$tau = 1/S^2
  }
  var_type = flashr:::handle_var_type(var_type, data)
  init_fn = flashr:::handle_init_fn(init_fn)
  ebnm_fn = flashr:::handle_ebnm_fn(ebnm_fn)
  ebnm_param = flashr:::handle_ebnm_param(ebnm_param, ebnm_fn, Kmax)

  if (verbose) {
    verbose_output = "odn" # objective, obj diff, nullcheck
  } else {
    verbose_output = ""
  }
  verbose_output = unlist(strsplit(verbose_output, split=NULL))
  stopping_rule = 'objective'

  history = list()
  prev_K = flashr:::flash_get_k(f)

  for(k in 1:Kmax){

    if (length(verbose_output) > 0) {
      flashr:::verbose_greedy_next_fl(prev_K + k, stopping_rule, tol)
    }


    old_f = f
    res = smooth_flash_r1(data,
                          f,
                          var_type,
                          init_fn,
                          tol,
                          ebnm_fn$l,
                          ebnm_param$l[[k]],
                          ebnm_fn$f,
                          ebnm_param$f[[k]],
                          filter.number,
                          family,
                          nullcheck,
                          maxiter,
                          verbose_output,
                          type)

    f = res$f
    next_history = res$history

    # Test whether the factor/loading combination is effectively zero.
    if (is_tiny_fl(f, flashr:::flash_get_k(f))) {
      # Remove zero factor as long as it doesn't create an empty object.
      if (flashr:::flash_get_k(f) > 1) {
        f = old_f
        next_history$zeroed_out = prev_K + k
      }
      history = c(history, list(next_history))
      break
    }

    history = c(history, list(next_history))

  }

  flash_object = flashr:::construct_flash_object(data = data,
                                        fit = f,
                                        history = history,
                                        f_init = f_init)

  return(flash_object)

}


smooth_flash_r1 = function(data,
                           f_init,
                           var_type,
                           init_fn,
                           tol,
                           ebnm_fn_l,
                           ebnm_param_l,
                           ebnm_fn_f,
                           ebnm_param_f,
                           filter.number,
                           family,
                           nullcheck,
                           maxiter,
                           verbose_output,
                           type
) {



  f = flashr:::add_factors_from_data(data, K = 1, f_init, init_fn)

  #browser()

  opt_res = smooth_flash_optimize_single_fl(data,
                                            f,
                                            flashr:::flash_get_k(f),
                                            var_type,
                                            tol,
                                            ebnm_fn_l,
                                            ebnm_param_l,
                                            ebnm_fn_f,
                                            ebnm_param_f,
                                            filter.number,
                                            family,
                                            maxiter,
                                            verbose_output,
                                            type
  )

  f = opt_res$f

  if (nullcheck) {
    null_res = flashr:::perform_nullcheck(data,
                                 f,
                                 flashr:::flash_get_k(f),
                                 var_type,
                                 verbose = ("n" %in% verbose_output))
    f = null_res$f
    # zeroed_out field is handled in flash_greedy_workhorse
  }

  return(list(f = f, history = opt_res$history))
}



smooth_flash_optimize_single_fl = function(data,f,k,
                                           var_type,
                                           tol,
                                           ebnm_fn_l,
                                           ebnm_param_l,
                                           ebnm_fn_f,
                                           ebnm_param_f,
                                           filter.number,
                                           family,
                                           maxiter,
                                           verbose_output,
                                           type){

  if (length(verbose_output) > 0) {
    flashr:::verbose_obj_table_header(verbose_output)
  }

  R2 = flashr:::flash_get_R2(data, f)
  # Expected residuals and squared residuals with factor k excluded:
  Rk = flashr:::flash_get_Rk(data, f, k)
  R2k = (R2 + 2 * outer(f$EL[, k], f$EF[, k]) * Rk
         - outer(f$EL2[, k], f$EF2[, k]))

  iter = 0
  diff = Inf
  diff_track = rep(NA, maxiter)
  obj_track = rep(NA, maxiter)

  while ((iter < maxiter) && (diff > tol)) {
    iter = iter + 1

    #if(iter==11){browser()}

    f = smooth_flash_update_single_lf(data,
                                      f,
                                      k,
                                      var_type,
                                      ebnm_fn_l,
                                      ebnm_param_l,
                                      ebnm_fn_f,
                                      ebnm_param_f,
                                      filter.number,
                                      family,
                                      Rk,
                                      R2,
                                      type)



    R2 = (R2k - 2 * outer(f$EL[, k], f$EF[, k]) * Rk
          + outer(f$EL2[, k], f$EF2[, k]))


    obj_track[iter] = (sum(unlist(f$KL_l)) + sum(unlist(f$KL_f)) +
                         flashr:::e_loglik_from_R2_and_tau(R2, f$tau, data))
    obj_diff = flashr:::calc_obj_diff(obj_track, iter)

    diff = flashr:::calc_diff("objective", obj_diff)
    #diff = abs(calc_diff("objective", obj_diff))
    diff_track[iter] = diff

    if (length(verbose_output) > 0) {
      flashr:::verbose_obj_table_entry(verbose_output,
                              iter,
                              obj_track[iter],
                              obj_diff,
                              max_chg_l=NULL,
                              max_chg_f=NULL,
                              f$gl[k],
                              f$gf[k])
    }

  }

  history = list(type = "greedy",
                 kset = k,
                 niter = iter,
                 obj_track = obj_track[1:iter],
                 diff_track = diff_track[1:iter])

  return(list(f = f, history = history))

}


smooth_flash_update_single_lf = function(data,f,k,
                                         var_type,
                                         ebnm_fn_l,
                                         ebnm_param_l,
                                         ebnm_fn_f,
                                         ebnm_param_f,
                                         filter.number,
                                         family,
                                         Rk = NULL,
                                         R2= NULL,
                                         type){

  #browser()

  if (is.null(R2)) {
    R2 = flashr:::flash_get_R2(data, f)
  }
  f$tau = flashr:::compute_precision(R2, data, var_type)

  if (is.null(Rk)) {
    Rk = flashr:::flash_get_Rk(data, f, k)
  }

  f = flashr:::flash_update_single_loading(data,
                                  f,
                                  k,
                                  ebnm_fn_l,
                                  ebnm_param_l,
                                  Rk,
                                  calc_obj = TRUE)

  f = smooth_flash_update_single_factor(data,
                                        f,
                                        k,
                                        ebnm_fn_f,
                                        ebnm_param_f,
                                        filter.number,
                                        family,
                                        Rk,
                                        calc_obj = TRUE,
                                        type)




  return(f)

}



smooth_flash_update_single_factor = function(data,f,k,
                                             ebnm_fn,
                                             ebnm_param,
                                             filter.number,family,
                                             Rk,calc_obj = TRUE,
                                             type){

  #browser()

  subset = which(!f$fixf[, k])
  #print(subset)
  any_fixed = any(f$fixf[, k])
  ebnm_args = flashr:::calc_ebnm_f_args(data, f, k, subset, any_fixed, Rk)
  if (is.null(ebnm_args)) {
    return(f)
  }

  if (!is.null(ebnm_param$warmstart)) {
    if (ebnm_param$warmstart) {
      if(length(f$gf) >= k){
        ebnm_param$g = f$gf[[k]]
      }
    }
    ebnm_param$warmstart = NULL
  }

  if(type=='wavelet'){
    a = smash_dwt(ebnm_args$x, ebnm_args$s,filter.number,family,ebnm_fn,ebnm_param)
    res = list(EX = a$mu.est,
               EX2 = a$mu.est.var+(a$mu.est)^2,
               g = NULL,
               KL = a$dKL)
  }
  if(type=='station'){
    a = smash.gaus(ebnm_args$x, ebnm_args$s,
                           filter.number=filter.number,
                           family=family,
                           return.loglr = TRUE,
                           post.var = TRUE)
    res = list(EX = a$mu.est,
               EX2 = a$mu.est.var+(a$mu.est)^2,
               g = NULL)
    # the obj in this case is not exact. So ELBO might not be increasing.
    if (calc_obj) {
      KL = a$loglik - flashr:::NM_posterior_e_loglik(ebnm_args$x, ebnm_args$s,
                                            res$EX, res$EX2)
      res$KL =KL
    }
  }

  #a = do.call(ebnm_fn, list(ebnm_args$x, ebnm_args$s, ebnm_param))




  if (!is.null(res)) {

    f$EF[subset, k] = res$EX
    f$EF2[subset, k] = res$EX2
    f$gf[[k]] = res$g
    if (calc_obj) {
      f$KL_f[[k]] = res$KL
    }

  }

  return(f)

}


