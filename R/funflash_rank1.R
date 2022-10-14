
#'@description  fit a rank 1 funflash model
#'@param Rk data matrix, or Rk, the residual matrix
#'@param res current fit of the model
#'@param k current k to be added
#'@return res
#'@export
funflash_rank1 = function(Rk,
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
                          type,
                          Y,
                          nlevel,
                          S2.type,
                          data){

  #s = res$s
  # initialize l and f
  n = nrow(Rk$Y)
  p = ncol(Rk$Y)
  tempY = get_Rk_data(data$orig.Y,res$Fe,res$a,res$b,res$EL,res$EF,k-1,data$filter.number,data$family,type,data$orig.p)
  init = init_fl_data(tempY,init_fn,data$filter.number,data$family,type)
  El = init$u*sqrt(init$d[1])
  Ef = init$v*sqrt(init$d[1])
  fe = init$e*sqrt(init$d[1])
  El2 = El^2
  Ef2 = Ef^2
  gl = list()
  gf = list()
  KL.l = 0
  KL.f = 0
  res = update_res(res,k,El,El2,Ef,Ef2,gl,gf,KL.l,KL.l,fe=fe)
  sigma2 = res$sigma2

  obj = -Inf
  for(i in 1:maxiter){

    # update sigma2
    if(var.type!='zero'){
      sigma2 = update_sigma2(Y,a,b,S2,res,var.type,S2.type,data$energy,data$energy.S2)
    }


    # update l

    ## formulate the x and s for ebnm

    #sl = (rowSums((a^2*tcrossprod(rep(1,n),b^2*Ef2)/(sigma2+S2))))^(-0.5)
    #xl = sl^2*rowSums(a*t(b*Ef*t(Rk))/(sigma2+S2))
    sl = (rowSums(a^2*tcrossprod(rep(1,n),c(Ef2,fe^2))/(sigma2+cbind(S2,data$energy.S2))))^(-0.5)
    xl = sl^2*rowSums(a*t(c(Ef,fe)*t(cbind(Rk$Y,Rk$energy)))/(sigma2+cbind(S2,data$energy.S2)))
    fit.l = do.call(ebnm_fn, list(xl, sl, ebnm_param))

    El = fit.l$postmean
    El2 = fit.l$postmean2
    gl = fit.l$fitted_g

    KL.l = fit.l$penloglik - NM_posterior_e_loglik(xl, sl,El, El2)

    if(sum(El2)==0){
      res = update_res(res,k,El,El2,Ef,Ef2,gl,gf,KL.l,KL.f,sigma2,fe)
      print('loading zeroed out')
      break
    }

    # update f

    # update wavelet coefficients
    KL.f = 0
    for(t in (0):(nlevel-1)){
      idx = get_wave_index(t,nlevel,type)
      yt = Rk$Y[,idx,drop=FALSE]
      #sf = (b[idx]^2*colSums(tcrossprod(a^2*El2,rep(1,length(idx)))/(sigma2[,idx,drop=FALSE]+S2[,idx,drop=FALSE])))^(-0.5)
      #xf = sf^2*b[idx]*colSums(yt*tcrossprod(a*El,rep(1,length(idx)))/(sigma2[,idx,drop=FALSE]+S2[,idx,drop=FALSE]))
      sf = (colSums(tcrossprod(a^2*El2,rep(1,length(idx)))/(sigma2[,idx,drop=FALSE]+S2[,idx,drop=FALSE])))^(-0.5)
      xf = sf^2*colSums(yt*tcrossprod(a*El,rep(1,length(idx)))/(sigma2[,idx,drop=FALSE]+S2[,idx,drop=FALSE]))

      fit.f = do.call(ebnm_fn, list(xf, sf, ebnm_param))
      Ef[idx] = fit.f$postmean
      Ef2[idx] = fit.f$postmean2
      gf[[t+1]] = fit.f$fitted_g
      KL.f = KL.f + fit.f$penloglik - NM_posterior_e_loglik(xf, sf,fit.f$postmean, fit.f$postmean2)
    }

    # update the fe
    fe = sum(Rk$energy*El/rowMeans(sigma2))/sum(El2/rowMeans(sigma2))

    if(sum(Ef2)==0){
      res = update_res(res,k,El,El2,Ef,Ef2,gl,gf,KL.l,KL.f,sigma2,fe)
      print('factor zeroed out')
      break
    }

    res = update_res(res,k,El,El2,Ef,Ef2,gl,gf,KL.l,KL.f,sigma2,fe)
    # evaluate objective function
    obj[i+1] = calc_objective_e(Y,a,b,res,S2,data$energy)
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
