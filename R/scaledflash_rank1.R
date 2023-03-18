
#'Fit rank 1 scaled flash
#'@export
scaledflash_rank1 = function(Rk,
                          a,
                          b,
                          S2,
                          res,
                          k,
                          var.type,
                          ebnm_fn,
                          ebnm_param,
                          maxiter=100,
                          tol=0.01,
                          verbose,
                          init_fn,
                          Y,
                          S2.type){

  #s = res$s
  # initialize l and f
  n = nrow(Rk)
  p = ncol(Rk)
  init = init_fl(Rk,init_fn)
  El = init$u*sqrt(init$d[1])/a
  Ef = init$v*sqrt(init$d[1])/b
  El2 = El^2
  Ef2 = Ef^2
  gl = list()
  gf = list()
  KL.l = 0
  KL.f = 0
  res = update_res(res,k,El,El2,Ef,Ef2,gl,gf,KL.l,KL.l,fe=0)
  sigma2 = res$sigma2

  obj = -Inf
  for(i in 1:maxiter){

    # update sigma2
    if(var.type!='zero'){
      sigma2 = update_sigma2(Y,a,b,S2,res,var.type,S2.type)
    }


    # update l

    ## formulate the x and s for ebnm

    #xl = rowSums(rep(1,n)%*%t(Ef) * Y)/sum(Ef2)
    sl = (rowSums((a^2*tcrossprod(rep(1,n),b^2*Ef2)/(sigma2+S2))))^(-0.5)
    xl = sl^2*rowSums(a*t(b*Ef*t(Rk))/(sigma2+S2))
    fit.l = do.call(ebnm_fn, list(xl, sl, ebnm_param))

    El = fit.l$postmean
    El2 = fit.l$postmean2
    gl = fit.l$fitted_g

    KL.l = fit.l$penloglik - NM_posterior_e_loglik(xl, sl,El, El2)

    if(sum(El2)==0){
      res = update_res(res,k,El,El2,Ef,Ef2,gl,gf,KL.l,KL.f,sigma2,0)
      print('loading zeroed out')
      break
    }

    # update f
    sf = (b^2*colSums(tcrossprod(a^2*El2,rep(1,p))/(sigma2+S2)))^(-0.5)
    xf = sf^2*b*colSums(Rk*tcrossprod(a*El,rep(1,p))/(sigma2+S2))
    fit.f = do.call(ebnm_fn, list(xf, sf, ebnm_param))
    Ef = fit.f$postmean
    Ef2 = fit.f$postmean2
    gf = fit.f$fitted_g
    KL.f = fit.f$penloglik - NM_posterior_e_loglik(xf, sf,Ef, Ef2)

    if(sum(Ef2)==0){
      res = update_res(res,k,El,El2,Ef,Ef2,gl,gf,KL.l,KL.f,sigma2,0)
      print('factor zeroed out')
      break
    }

    res = update_res(res,k,El,El2,Ef,Ef2,gl,gf,KL.l,KL.f,sigma2,0)
    # evaluate objective function
    obj[i+1] = calc_objective(Y,a,b,res,S2)
    if(verbose){
      print(paste('Iteration ', i, ': obj ', round(obj[i+1],3)))
    }

    if((obj[i+1]-obj[i])<0){
      print('An iteration decreased the objective')
    }

    if((obj[i+1]-obj[i])<=tol){
      break
    }


  }

  return(res)


}
