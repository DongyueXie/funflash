#'@title perform biwhitening of a matrix
#'@param Y A Poisson matrix
#'@export
biwhitening = function(Y){
  out = Sinkhorn_Knopp(Y)
  u = sqrt(out$x)
  v = sqrt(out$y)
  # Y_tilde = sqrt(u) * Y%*%diag(v)
  Y_tilde = u * Y%*%diag(v)
  return(list(Y = Y_tilde,u=u,v=v))
}


#'@title Sinkhorn_Knopp for matrix row and col scaling
#'@param A nonnegative matrix
#'@param rs target row sums, a vector of length m
#'@param cs target col sums, a vector of length n
#'@param tol tol to stop the iterations
#'@export
Sinkhorn_Knopp = function(A,rs=NULL,cs=NULL,tol=1e-5,maxiter=100){

  m = nrow(A)
  n = ncol(A)

  if(is.null(rs)){
    rs = rep(n,m)
  }
  if(is.null(cs)){
    cs = rep(m,n)
  }

  x = rep(1,m)
  y = rep(1,n)


  niter = 0
  while(niter<=maxiter) {
    Ax = c(crossprod(A,x))
    if(max(abs(y*Ax))<= tol){
      break
    }
    y = cs/Ax
    Ay = c(A%*%y)
    x = rs/Ay
    #if(max(abs(x*Ay-rs))<=tol & max(abs(y*c(crossprod(A,x))-cs))<= tol){
    #  break
    #}
    niter = niter + 1
  }

  return(list(x=as.vector(x),y=as.vector(y)))

}
