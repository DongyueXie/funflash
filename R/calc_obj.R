# @title Expected log likelihood for normal means model.
#
# @description The likelihood is for x | theta ~ N(theta, s^2);
#   The expectation is taken over the posterior on theta.
#
# @param x Observations in normal means.
#
# @param s Standard errors of x.
#
# @param Et Posterior mean of theta.
#
# @param Et2 Posterior second moment of theta.
#
NM_posterior_e_loglik = function(x, s, Et, Et2) {
  # Deal with infinite SEs:
  idx = is.finite(s)
  x = x[idx]
  s = s[idx]
  Et = Et[idx]
  Et2 = Et2[idx]
  return(-0.5 * sum(log(2*pi*s^2) + (1/s^2) * (Et2 - 2*x*Et + x^2)))
}


#'@description evaluate objective function
calc_objective = function(Y,a,b,res,S2){
  # is(is.null(S2)){
  #   obj = Eloglik(Y,a,b,res$sigma2,res$EL,res$EF,res$EL2,res$EF2) + sum(res$KL_L) + sum(res$KL_F)
  # }else{
  #   obj = Eloglik(Y,a,b,res$sigma2+S2,res$EL,res$EF,res$EL2,res$EF2) + sum(res$KL_L) + sum(res$KL_F)
  # }
  obj = Eloglik(Y,a,b,res$sigma2+S2,res$EL,res$EF,res$EL2,res$EF2) + sum(res$KL_L) + sum(res$KL_F)
  obj
}

#'@description This function calculates the expected log likelihood wrt q(L,F).
#'@param Y the data matrix
#'@param Sigma a matrix of variance of random errors. [sigma^2 sigma^2 ... sigma^2]
Eloglik = function(Y,a,b,Sigma,EL,EF,EL2,EF2){
  -0.5 * sum(log(2 * pi * Sigma) + 1/Sigma * get_R2(Y,a,b,EL,EF,EL2,EF2))
}



#'@description evaluate objective function
calc_objective_e = function(Y,a,b,res,S2,energy){
  obj = Eloglik_e(Y,a,b,res$sigma2[,-ncol(res$sigma2)]+S2,res$EL,res$EF,res$EL2,res$EF2,energy,res$Fe) + sum(res$KL_L) + sum(res$KL_F)
  obj
}

#'@description This function calculates the expected log likelihood wrt q(L,F).
#'@param Y the data matrix
#'@param Sigma a matrix of variance of random errors. [sigma^2 sigma^2 ... sigma^2]
Eloglik_e = function(Y,a,b,Sigma,EL,EF,EL2,EF2,energy,Fe){
  R2 = get_R2_e(Y,energy,Fe,a,b,EL,EF,EL2,EF2)
  -0.5 * sum(log(2 * pi * Sigma) + 1/Sigma * R2$T2_all)
}

