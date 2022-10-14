#'Check if the added factor is null
null_check = function(res,k,verbose,Y,a,b,S2,var.type,S2.type){
  # delete the kth factor
  res0 = delete_factor(res,k,Y,a,b,S2,var.type,S2.type)
  # compare the objective function
  obj0 = calc_objective(Y,a,b,res0,S2)
  obj1 = calc_objective(Y,a,b,res,S2)

  if(obj0>obj1){
    res = res0
  }
  if(verbose){
    if(obj1>obj0){
      print(paste('Deleting factor ',k, ' decreases objective by ', round(obj1-obj0,3)))
    }else{
      print(paste('Deleting factor ',k, ' increases objective by ', round(obj0-obj1,3)))
    }
  }
  res
}

#'Delete the kth dimension from res object
delete_factor = function(res,k,Y,a,b,S2,var.type,S2.type){
  res$EL[,k] = 0
  res$EF[,k] = 0
  res$EL2[,k] = 0
  res$EF2[,k] = 0
  res$gL[[k]] = list(NULL)
  res$gF[[k]] = list(NULL)
  res$KL_L[k] = 0
  res$KL_F[k] = 0
  if(var.type!='zero'){
    res$sigma2 = update_sigma2(Y,a,b,S2,res,var.type,S2.type)
  }
  res
}

#'Check if the added loading and factors are very small. If so, break the iteration.
is_tiny_fl = function (f, k, tol = 1e-08) {
  return(sum(f$EL[, k]^2) * sum(f$EF[, k]^2) < tol)
}
