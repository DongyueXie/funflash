
#'Initialize a rank 1 decomposition
init_fl = function(Y,init_fn){
  if(init_fn=='udv_si'){
    out = udv_si(Y,1)
  }
  if(init_fn=='udv_si_svd'){
    out = udv_si_svd(Y,1)
  }
  if(init_fn=='udv_svd'){
    out = udv_svd(Y,1)
  }
  if(init_fn=='udv_random'){
    out = udv_random(Y,1)
  }
  if(init_fn == 'nnmf_r1'){
    out = nnmf_r1(Y,1)
  }

  #out = do.call(init_fn,list(Y,1))
  out
}

nnmf_r1 = function(Y,K = 1){
  res = NNLM::nnmf(Y,K,loss='mse',verbose = 0)
  u = as.vector(res$W)
  v = as.vector(res$H)
  d = sum(v) * sum(u)
  u = u/sum(u)
  v = v/sum(v)
  return(list(u=u,v=v,d=d))
}

init_fl_data = function(Y,init_fn,filter.number,family,type){
  if(init_fn=='udv_si'){
    out = udv_si(Y,1)
  }
  if(init_fn=='udv_si_svd'){
    out = udv_si_svd(Y,1)
  }
  if(init_fn=='udv_svd'){
    out = udv_svd(Y,1)
  }
  if(init_fn=='udv_random'){
    out = udv_random(Y,1)
  }
  temp = wd(out$v,filter.number=filter.number,family=family,type=type)
  out$e = sum(out$v)/sqrt(length(out$v))
  out$v = temp$D
  return(out)
}

# @title udv_si
#
# @description Provides a simple wrapper to \code{softImpute} to
#   provide a rank 1 initialization. Uses \code{type = "als"} option.
#
# @param Y An n by p matrix.
#
# @param K Number of factors to use.
#
# @return A list with components (u,d,v).
#
#' @importFrom softImpute softImpute
#'
udv_si = function(Y, K = 1) {
  suppressWarnings(
    res <- softImpute(Y, rank.max = K, type = "als", lambda = 0)
  )
  return(res)
}


# @title udv_si_svd
#
# @description provides a simple wrapper to \code{softImpute} to
#   provide a rank 1 initialization. Uses \code{type = "svd"} option.
#
# @inherit udv_si
#
#' @importFrom softImpute softImpute
#'
udv_si_svd = function(Y, K = 1) {
  suppressWarnings(
    res <- softImpute(Y, rank.max = K, type = "svd", lambda = 0)
  )
  return(res)
}


# @title udv_svd
#
# @description Provides a simple wrapper to svd.
#
# @inherit udv_si
#
udv_svd = function (Y, K = 1) {
  svd(Y, K, K)
}


# @title udv_random
#
# @description Provides a random initialization of factors.
#
# @inherit udv_si
#
# @return A list with components (u,d,v), with elements of u and v
#   i.i.d. N(0,1).
#
#' @importFrom stats rnorm
#'
udv_random = function (Y, K = 1) {
  n = nrow(Y)
  p = ncol(Y)
  return(list(u = matrix(rnorm(n * K), ncol = K),
              d = 1,
              v = matrix(rnorm(p * K), ncol = K)))
}


