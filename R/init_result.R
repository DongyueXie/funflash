
#'@description initialize the output list
init_res = function(n,p,Kmax,sigma2){
  if(is.null(sigma2)){
    sigma2 = matrix(0,nrow=n,ncol=p+1)
  }else{
    if(length(sigma2)==1){
      sigma2 = matrix(sigma2,nrow=n,ncol=p)
    }else if(length(sigma2)==n){
      sigma2 = tcrossprod(sigma2,rep(1,p))
    }else if(length(sigma2==p)){
      sigma2 = tcrossprod(rep(1,n),sigma2)
    }
  }
  return(list(EL = matrix(0,nrow=n,ncol=Kmax),
              EF = matrix(0,nrow=p,ncol=Kmax),
              EL2 = matrix(0,nrow=n,ncol=Kmax),
              EF2 = matrix(0,nrow=p,ncol=Kmax),
              Fe = rep(0,Kmax),
              gL = list(),
              gF = list(),
              KL_L = rep(0,Kmax),
              KL_F = rep(0,Kmax),
              sigma2 = sigma2))
}

#'@description  update the output list
update_res = function(res,k,El,El2,Ef,Ef2,
                      gl,gf,
                      KL.l,KL.f,
                      sigma2=NULL,
                      fe){
  res$EL[,k] = El
  res$EL2[,k] = El2
  res$EF[,k] = Ef
  res$EF2[,k] = Ef2
  res$gL[[k]] = gl
  res$gF[[k]] = gf
  res$KL_L[k] = KL.l
  res$KL_F[k] = KL.f
  res$Fe[k] = fe
  if(!is.null(sigma2)){
    res$sigma2 = sigma2
  }
  res
}
