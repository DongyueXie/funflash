
#' Calculate R2 for updating the variance
get_R2 = function(Y,a,b,EL,EF,EL2,EF2){
  (Y-tcrossprod(a*EL,b*EF))^2+tcrossprod(a^2*EL2,b^2*EF2) - tcrossprod(a^2*EL^2,b^2*EF^2)
}

#' Calculate R2 for updating the variance
get_R2_e = function(Y,energy,Fe,a,b,EL,EF,EL2,EF2){
  R2 = (Y-tcrossprod(a*EL,b*EF))^2+tcrossprod(a^2*EL2,b^2*EF2) - tcrossprod(a^2*EL^2,b^2*EF^2)
  R2_e = (energy - a*EL%*%(Fe))^2 + a^2*EL2%*%(Fe^2) - a^2*EL^2%*%(Fe^2)
  return(list(R2=R2,R2_e=R2_e,R2_all=cbind(R2,R2_e)))
}

#'@description get the residual for model fitting, kth factor and loadings
get_Rk = function(Y,energy,Fe,a,b,EL,EF){
  if(is.null(EL)&is.null(EF)){
    return(list(Y=Y,energy=energy))
  }else{
    return(list(Y=Y - tcrossprod(a*EL,b*EF),energy = energy - a*EL%*%(Fe)))
  }
}

#'@description get the residual for model fitting, kth factor and loadings
get_Rk0 = function(Y,a,b,EL,EF){
  if(is.null(EL)&is.null(EF)){
    return(Y)
  }else{
    return(Y - tcrossprod(a*EL,b*EF))
  }
}

#'@description get the residual for model fitting, kth factor and loadings
#'@param Y original Y in the data space
#'@param p the ncol of original data Y
#'@param k number of factors added so far
get_Rk_data = function(Y,Fe,a,b,EL,EF,k,filter.number,family,type,p){
  if(is.null(EL)&is.null(EF)|k==0){
    return(Y)
  }else{
    EF = inv_dwt_wc(EF,Fe,p,filter.number,family,type)
    return(Y - tcrossprod(EL,EF))
  }
}

