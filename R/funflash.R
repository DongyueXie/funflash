#'@title Fit smoothing flash
#'@param data a data obj from funflash_set_data()
#'@param var.type "by_column","by_row", "constant", "zero".
#'@param Kmax max number of factors
#'@param tol tol for stopping the iterations
#'@param maxiter max number of iteration
#'@param ebnm_fn 'ebnm_pn', 'ebnm_pl', 'ebnm_ash'
#'@param init_fn 'udv_si',..., see ?flashr
#'@param nullcheck whether perform null check
#'@param sigma2 if not null and var.type='zero', treat it as given. Can be used as test purpose.
#'@return a list of
#'  \item{ldf}{a list of d, El, Ef}
#'  \item{Fb}{E(Fb) before dividing b}
#'  \item{nfactors}{number of factors}
#'  \item{pve}{proportion/percentage variance explained}
#'  \item{fitted.values}{aELEF'b}
#'  \item{fit}{result list, before transforming back to data space}
#'  \item{objective}{objective funtion value}
#'@export
funflash = function(data,
                    var.type = "by_row",
                    Kmax=10,
                    tol=0.01,
                    maxiter = 1000,
                    ebnm_fn = 'ebnm_pn',
                    init_fn = 'udv_si',
                    ebnm_param=NULL,
                    verbose=TRUE,
                    nullcheck=TRUE,
                    sigma2 = NULL,
                    seed=12345){

  set.seed(seed)

  n = nrow(data$Y)
  p = ncol(data$Y)

  res = init_res(n,p,Kmax,sigma2)
  #res$energy = data$energy

  for(k in 1:Kmax){
    res_old = res
    Rk = get_Rk(data$Y,data$energy,res$Fe,data$a,data$b,res$EL,res$EF)
    if(verbose){
      print(paste("Fitting dimension ", k))
    }
    res = funflash_rank1(Rk,
                         a = data$a,
                         b = data$b,
                         data$S2,
                         res,
                         k,
                         var.type,
                         ebnm_fn=ebnm_fn,
                         ebnm_param=ebnm_param,
                         maxiter=maxiter,
                         tol=tol,
                         verbose=verbose,
                         init_fn=init_fn,
                         type=data$type,
                         Y=data$Y,
                         nlevel=data$nlevel,
                         S2.type = data$S2.type,
                         data=data)
    #est_var=est_var)
    if(nullcheck){
      if(verbose){
        print("Performing nullcheck")
      }
      res = null_check(res,k,verbose,data$Y,data$a,data$b,data$S2,var.type,data$S2.type)
    }

    if(is_tiny_fl(res,k)){
      if(k>1){
        res=res_old
      }
      break
    }
    # check to stop
  }

  # transform back to data space
  # also deal with non-power 2: transform it back
  #browser()
  out = construct_object_e(res,data,var.type)
  out

}





