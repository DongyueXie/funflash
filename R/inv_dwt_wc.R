#'@description get the NDWT wavelet coefficients
#'@param p the dimension of original Y, possible after reflect.
#'@param Y input data matrix
inv_dwt_wc = function(EF,f_e,p,filter.number,family,type){
  K = ncol(EF)
  FF = matrix(nrow=p,ncol=K)
  if(type=='wavelet'){

    for(k in 1:K){
      temp0 = wavethresh::wd(rep(0,p),filter.number=filter.number,family=family,type='wavelet')
      temp0$D = EF[,k]
      temp0 = putC(temp0,0,f_e[k])
      FF[,k] = wr(temp0)
    }

  }else if(type == 'station'){

    for(k in 1:K){
      temp0 = wavethresh::wd(rep(0,p),filter.number=filter.number,family=family,type='station')
      temp0$D = EF[,k]
      temp0 = putC(temp0,0,rep(f_e[k],p))
      FF[,k] = AvBasis(convert(temp0))
    }

  }else{
    stop('type should be either wavelet or station')
  }
  return(FF)
}
