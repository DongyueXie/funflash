#'@title Fit scaled flash
#'@param Y a data obj from funflash_set_data()
#'@param a row scaling factors
#'@param b column scaling factors
#'@param S2 Known variances. Can be NULL, a scalar, a vector or a matrix.
#'@param var.type whether assume F is smooth.
#'@export
scaledflash = function(Y,
                    a=NULL,
                    b=NULL,
                    S2 = NULL,
                    var.type = 'by_row',
                    Kmax=50,
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
  n = nrow(Y)
  p = ncol(Y)

  if(is.null(a)){
    a = rep(1,n)
  }
  if(is.null(b)){
    b = rep(1,p)
  }

  if(is.null(S2)){
    S2 = matrix(0,nrow=n,ncol=p)
    S2.type='zero'
  }else{
    if(length(S2)==n){
      S2 = tcrossprod(S2,rep(1,n))
    }else if(length(S2)==p){
      S2 = tcrossprod(rep(1,n),S2)
    }else if(length(S2)==1){
      S2 = matrix(S2,nrow=n,ncol=p)
    }
    S2.type = 'matrix'
  }



  res = init_res(n,p,Kmax,sigma2)
  #res$energy = data$energy

  for(k in 1:Kmax){
    res_old = res
    Rk = get_Rk(Y,a,b,res$EL,res$EF)
    if(verbose){
      print(paste("Fitting dimension ", k))
    }
    res = scaledflash_rank1(Rk,
                         a = a,
                         b = b,
                         S2,
                         res,
                         k,
                         var.type,
                         ebnm_fn=ebnm_fn,
                         ebnm_param=ebnm_param,
                         maxiter=maxiter,
                         tol=tol,
                         verbose=verbose,
                         init_fn=init_fn,
                         Y=Y,
                         S2.type = S2.type)
    #est_var=est_var)
    if(nullcheck){
      if(verbose){
        print("Performing nullcheck")
      }
      res = null_check(res,k,verbose,Y,a,b,S2,var.type,S2.type)
    }

    if(is_tiny_fl(res,k)){
      if(k>1){
        res=res_old
      }
      break
    }
    # check to stop
  }


  # remove factor and loadings are zero
  rm.idx = which(colSums(res$EL2)==0)
  if(length(rm.idx)>0){
    res$EL = res$EL[,-rm.idx,drop=FALSE]
    res$EF = res$EF[,-rm.idx,drop=FALSE]
    res$EL2 = res$EL2[,-rm.idx,drop=FALSE]
    res$EF2 = res$EF2[,-rm.idx,drop=FALSE]
    res$KL_L = res$KL_L[-rm.idx]
    res$KL_F = res$KL_F[-rm.idx]
  }

  d = sqrt(colSums(res$EL^2) * colSums(res$EF^2))
  fitted_values = tcrossprod(res$EL,res$EF)

  FF = scale(res$EF,scale = sqrt(colSums(res$EF^2)),center=FALSE)
  L = scale(res$EL,scale = sqrt(colSums(res$EL^2)),center=FALSE)

  ldf = list(d=d,l=L,f=FF)

  nfactors = length(d)
  pve = d^2/(sum(d^2)+sum(S2+res$sigma2))
  objective = calc_objective(Y,a,b,res,S2)



  return(list(ldf = ldf,
              nfactors = nfactors,
              pve = pve,
              fitted_values = fitted_values,
              fit = res,
              objective = objective))



}





