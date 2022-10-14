#'@description construct object for returning results.
construct_object = function(res,data,var.type){
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

  # take the energy into account
  # here E(energy|EL) = EL*f_e; var(energy|EL) = diag(s^2)
  # perform a weighted least square to estimate f_e
  #browser()
  if(var.type=='by_column'){
    f_e = lm(y~.+0,data.frame(y=data$energy,X = data$a*res$EL))
    # in this case, we first estimate f_e by ols, then estimate the sigma2 for energy.
    energy.sigma2 = optimize_sigma2(data$energy^2,data$energy.S2)
  }else{
    energy.sigma2 = res$sigma2[,1]
    f_e = lm(y~.+0,data.frame(y=data$energy,X = data$a*res$EL),weights = 1/(data$energy.S2+energy.sigma2))
  }
  f_e = f_e$coefficients

  # calculate pve in wavelet space.
  if(data$type=='wavelet'){
    d = sqrt(colSums(data$a^2*res$EL^2) * (colSums(res$EF^2)+f_e^2))
    pve = d^2/(sum(d^2)+sum(energy.sigma2+data$energy.S2+data$S2+res$sigma2))
  }else if(data$type == 'station'){
    # calculate colSums
    p.temp = 2^data$nlevel
    d = sqrt(colSums(data$a^2*res$EL^2) * (colSums(res$EF^2)+f_e^2*p.temp))
    #pve = d^2/(sum(d^2)+sum(energy.sigma2*p.temp+data$energy.S2*p.temp+data$S2+res$sigma2))
    pve = 'to be done - need revise'
  }


  # transform F back to data space
  FF = inv_dwt_wc(res$EF,f_e,data$orig.p,data$filter.number,data$family,data$type)
  FF = FF[data$orig.idx,]

  # divide FF by the original b to get true F
  Fb = FF
  FF = FF/data$orig.b

  #d = sqrt(colSums(res$EL^2) * colSums(FF^2))
  fitted_values = tcrossprod(data$a*res$EL,data$orig.b*FF)

  FF = scale(FF,scale = sqrt(colSums(FF^2)),center=FALSE)
  L = scale(res$EL,scale = sqrt(colSums(res$EL^2)),center=FALSE)

  ldf = list(d=d,l=L,f=FF)

  nfactors = length(d)
  objective = calc_objective(data$Y,data$a,data$b,res,data$S2)



  return(list(ldf = ldf,
              Fb = Fb,
              nfactors = nfactors,
              pve = pve,
              fitted_values = fitted_values,
              fit = res,
              objective = objective))

}



#'@description construct object for returning results.
construct_object_e = function(res,data,var.type){
  # remove factor and loadings are zero
  rm.idx = which(colSums(res$EL2)==0)
  if(length(rm.idx)>0){
    res$EL = res$EL[,-rm.idx,drop=FALSE]
    res$EF = res$EF[,-rm.idx,drop=FALSE]
    res$EL2 = res$EL2[,-rm.idx,drop=FALSE]
    res$EF2 = res$EF2[,-rm.idx,drop=FALSE]
    res$KL_L = res$KL_L[-rm.idx]
    res$KL_F = res$KL_F[-rm.idx]
    res$Fe = res$Fe[-rm.idx]
  }


  # calculate pve in wavelet space.
  if(data$type=='wavelet'){
    d = sqrt(colSums(data$a^2*res$EL^2) * (colSums(res$EF^2)+res$Fe^2))
    #pve = d^2/(sum(d^2)+sum(energy.sigma2+data$energy.S2+data$S2+res$sigma2))
  }else if(data$type == 'station'){
    # calculate colSums
    p.temp = 2^data$nlevel
    d = sqrt(colSums(data$a^2*res$EL^2) * (colSums(res$EF^2)+res$Fe^2*p.temp))
    #pve = d^2/(sum(d^2)+sum(energy.sigma2*p.temp+data$energy.S2*p.temp+data$S2+res$sigma2))
    #pve = 'to be done - need revise'
  }


  # transform F back to data space
  FF = inv_dwt_wc(res$EF,res$Fe,data$orig.p,data$filter.number,data$family,data$type)
  FF = FF[data$orig.idx,]

  # divide FF by the original b to get true F
  #Fb = FF
  #FF = FF/data$orig.b

  #d = sqrt(colSums(res$EL^2) * colSums(FF^2))
  #fitted_values = tcrossprod(data$a*res$EL,data$orig.b*FF)

  FF = scale(FF,scale = sqrt(colSums(FF^2)),center=FALSE)
  L = scale(res$EL,scale = sqrt(colSums(res$EL^2)),center=FALSE)

  ldf = list(d=d,l=L,f=FF)

  nfactors = length(d)
  #objective = calc_objective(data$Y,data$a,data$b,res,data$S2)



  return(list(ldf = ldf,
              nfactors = nfactors,
              fit = res))

}
