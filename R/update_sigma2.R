
#'Update the variance
#'@param Y data matrix
#'@param a,b scaling factors
#'@return updated sigma2, a matrix
update_sigma2 = function(Y,a,b,S2,res,var.type,S2.type,energy,energy.S2){

  n = nrow(Y)
  p = ncol(Y)

  # todo: modify other var.type to take energy into account.
  if(var.type=='by_row'){
    if(S2.type=='zero'){
      sigma2 = c(rowMeans(get_R2(Y,a,b,res$EL,res$EF,res$EL2,res$EF2)))
    }else{
      R2 = get_R2(Y,a,b,res$EL,res$EF,res$EL2,res$EF2)
      #S2_e = cbind(S2,energy.S2)
      sigma2 = rep(0,n)
      for(r in 1:n){
        sigma2[r] = optimize_sigma2(R2[r,],S2[r,])
      }
    }
    sigma2 = tcrossprod(sigma2,rep(1,p))
  }else if(var.type=='by_column'){
    if(S2.type=='zero'){
      sigma2 = c(colMeans(get_R2(Y,a,b,res$EL,res$EF,res$EL2,res$EF2)))
    }else{
      R2 = get_R2(Y,a,b,res$EL,res$EF,res$EL2,res$EF2)
      sigma2 = rep(0,p)
      for(r in 1:p){
        sigma2[r] = optimize_sigma2(R2[,r],S2[,r])
      }
    }
    sigma2 = tcrossprod(rep(1,n),sigma2)
  }else if(var.type=='constant'){
    if(S2.type=='zero'){
      sigma2 = c(mean(get_R2(Y,a,b,res$EL,res$EF,res$EL2,res$EF2)))
    }else{
      R2 = get_R2(Y,a,b,res$EL,res$EF,res$EL2,res$EF2)
      sigma2 = optimize_sigma2(c(R2),c(S2))
    }
    sigma2 = matrix(sigma2,nrow=n,ncol=p)
  }

  return(sigma2)

}
