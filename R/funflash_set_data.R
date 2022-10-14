#'Set data for funflash
#'@param Y data matrix
#'@param a row scaling factors
#'@param b column scaling factors, should be smooth
#'@param S2 known variance. Can be a scalar, a vector or a matrix
#'@param family,filter.number wavelet family and filter number as in wavethresh package
#'@param type type of wavelet transformation, "wavelet" or "station". Choose 'station' could give better fit of the smooth factor but would introduce redundant factors.
#'@param reflect.data whether reflect the data matrix. Setting this to TRUE would give better fit to smooth F at boundary, but would likely introduce redundant factors.
#'@return a list of:
#'  \item{Y}{wavelet coefficients matrix, excluding the energy term}
#'  \item{a}{row scaling factor}
#'  \item{b}{a 1 vector. Because we need to DWT the bF together. Treat bF as the F and will divide b to get original F}
#'  \item{S2,S2.type}{S2 the known variance matrix, and its type}
#'  \item{orig.idx}{original index of Y before making it a power of 2 and possible reflection}
#'  \item{nlevel}{number of levels}
#'  \item{energy}{scaled summation of each row of Y}
#'  \item{orig.p}{number of column of Y after reflect() before DWT}
#'  \item{orig.b}{the input b}
#'  \item{orig.Y}{the original Y matrix(after reflection/power of 2)}
#'  \item{energy.S2}{the S2 for energy vector}
#'
#'@export
#'@import wavethresh
funflash_set_data = function(Y,
                             a = NULL,
                             b = NULL,
                             S2 = NULL,
                             family="DaubExPhase",
                             filter.number=1,
                             type='wavelet',
                             reflect.data = FALSE){

  input = list(Y=Y,
               a=a,
               b=b,
               S2=S2,
               family=family,
               filter.number=filter.number,
               type=type,
               reflect.data=reflect.data)
  # if S2 is a vector of length ncol(Y), then we get diag(WS2W') as the variance for wavelet coefficients
  # S2 is still a vector, finest level comes first.
  # If S2 is matrix, we do the transformation for each row of S2
  if(!is.null(S2)){
    if(length(S2) == ncol(Y)){
      S2 = reflect(S2,reflect.data)
      n = length(S2$x)
      if(type=='wavelet'){
        W2 = t(GenW(n = n,filter.number = filter.number,family = family)^2)
        S2 = (apply((rep(1, n) %o% S2$x) * W2, 1, sum))[-1]
      }else if(type=='station'){
        W2 = ndwt.mat(n = n,filter.number = filter.number,family = family)
        S2 = apply((rep(1, n * log2(n)) %o% S2$x) * W2, 1, sum)
      }
      energy.S2 = rep(mean(S2),nrow(Y))
      S2.type = 'matrix'
    }else if(length(S2) == prod(dim(Y))){
      S2 = reflect(S2,reflect.data)
      n = ncol(S2$x)
      if(type=='wavelet'){
        W2 = t(GenW(n = n,filter.number = filter.number,family = family)^2)
        S2 = t(apply(S2$x,1,function(z){(apply((rep(1, n) %o% z) * W2, 1, sum))[-1]}))
      }else if(type=='station'){
        W2 = ndwt.mat(n = n,filter.number = filter.number,family = family)
        S2 = t(apply(S2$x,1,function(z){apply((rep(1, n*log2(n)) %o% z) * W2, 1, sum)}))
      }
      energy.S2 = rowMeans(S2)
      S2.type = 'matrix'
    }else{
      stop('check dimension of S2')
    }
  }

  if(is.null(a)){
    a = rep(1,nrow(Y))
  }
  if(is.null(b)){
    orig.b = rep(1,ncol(Y))
  }else{
    orig.b = b
  }


  Y = reflect(Y,reflect.data)
  orig.idx = Y$idx
  Y = Y$x
  p = ncol(Y)
  nlevel = log2(p)

  # rows are wavelet coefficients. finest level comes first, the summation/energy term is not included.
  Yw = get_dwt_wc(Y,filter.number,family,type)

  if(length(S2)==1){
    energy.S2 = rep(S2,nrow(Yw))
    S2 = matrix(S2,nrow=nrow(Yw),ncol=ncol(Yw))
    S2.type = 'constant'
  }else if(length(S2)==nrow(Y)){
    energy.S2 = S2
    S2 = tcrossprod(S2,rep(1,ncol(Yw)))
    S2.type = 'by_row'
  }

  # ebnm cannot deal with ebnm(0,s) and it gives an error. So we set zeros to .Machine$double.eps. This happens when the coarsest level is 0 -- signal is reflected
  if(type=='wavelet'){
    if(any(Yw[,ncol(Yw)]==0)){
      Yw[,ncol(Yw)] = .Machine$double.eps
    }
  }


  # energy is the scaled summation of each row of Y
  energy = rowSums(Y)/sqrt(p)


  if(is.null(S2)){
    S2 = matrix(0,nrow=nrow(Yw),ncol=ncol(Yw))
    energy.S2 = rep(0,nrow(Yw))
    S2.type = 'zero'
  }

  ## ** note **##
  ## Here S2 is a matrix for every case
  ## Can be optimized later for lower memory usage.

  return(list(Y = Yw,
              a = a,
              b = rep(1,ncol(Yw)),
              S2 = S2,
              S2.type = S2.type,
              orig.idx=orig.idx,
              nlevel=nlevel,
              energy = energy,
              family=family,
              filter.number=filter.number,
              type=type,
              orig.p = p,
              orig.b = orig.b,
              orig.Y = Y,
              energy.S2=energy.S2,
              input=input))

}

# Modified wd() function from package 'wavethresh' to only return the
# detail coefficients. It returns the detail coefficients of a
# standard wavelet decomposition.
#
#' @importFrom wavethresh wd
wd.D = function (data, filter.number = 10, family = "DaubLeAsymm",
                 type = "wavelet", bc = "periodic", verbose = FALSE,
                 min.scale = 0, precond = TRUE) {
  l = wd(data = data, filter.number = filter.number,
                     family = family, type = type, bc = bc,
                     verbose = verbose, min.scale = min.scale,
                     precond = precond)
  return(l$D)
}

# @description Computes the non-decimated wavelet transform matrix for
#   a given basis.
# @param n The sample size. Must be a power of 2.
# @param filter.number Specifies the type of wavelet basis used.
# @param family Specifies the type of wavelet basis used.
# @return The NDWT matrix for the specified basis, with the entries squared.
ndwt.mat = function (n, filter.number, family) {
  J = log2(n)
  X = diag(rep(1, n))
  W = matrix(0, n * J, n)
  W = apply(X, 1, wd.D, filter.number = filter.number, family = family,
            type = "station")
  return(W^2)
}

#'@description get the DWT wavelet coefficients
#'@param Y input data matrix
#'@return wavelet coefficients. finest level comes first, the summation term is not included.
#'@import wavethresh
get_dwt_wc = function(Y,filter.number,family,type){
  Yw = apply(Y,1,function(z){
    w.d = wd(z,filter.number = filter.number,family = family,type=type)
    w.d$D
  })
  return(t(Yw))
}


#' @title Reflect and extend a vector, or each row of a matrix.
#'
#' @description Extends the vector to have length a power of 2 (if not already a power of 2) and then
#'   reflects it about its right end.
#'
#' @details The vector x is first reflected about both its left and right ends, by (roughly) the same
#' amount each end, to make its length a power of 2 (if the length of x is already a power of 2 this step is skipped).
#' Then the resulting vector is reflected about its right end to create a vector
#' that is both circular (left and right ends are the same) and a power of 2.
#'
#' @param x An n-vector.
#'
#' @return A list with two list elements: \code{"x"} containing the
#'   extended-and-reflected signal; and \code{"idx"} containing the indices of the
#'   original signal.
reflect <- function (x,reflect.data = TRUE) {

  if(is.matrix(x)){
    n = ncol(x)
    J = log2(n)
    if ((J%%1) == 0) {
      # if J is an integer, i.e. n is a power of 2.
      if(reflect.data){
        x = cbind(x, x[,n:1])
      }
      return(list(x=x, idx = 1:n))
    } else {
      n.ext = 2^ceiling(J)
      lnum = round((n.ext - n)/2)
      rnum = n.ext - n - lnum
      if (lnum == 0) {
        x.lmir = NULL
      } else {
        x.lmir = x[,lnum:1]
      }
      if (rnum == 0) {
        x.rmir = NULL
      } else {
        x.rmir = x[,n:(n - rnum + 1)]
      }
      x = cbind(x.lmir, x, x.rmir)
      if(reflect.data){
        x = cbind(x, x[,n.ext:1])
      }
      return(list(x = x, idx = (lnum + 1):(lnum + n)))
    }
  }else{
    n = length(x)
    J = log2(n)
    if ((J%%1) == 0) {

      # if J is an integer, i.e. n is a power of 2.
      if(reflect.data){
        x = c(x, x[n:1])
      }

      return(list(x=x, idx = 1:n))
    } else {
      n.ext = 2^ceiling(J)
      lnum = round((n.ext - n)/2)
      rnum = n.ext - n - lnum
      if (lnum == 0) {
        x.lmir = NULL
      } else {
        x.lmir = x[lnum:1]
      }
      if (rnum == 0) {
        x.rmir = NULL
      } else {
        x.rmir = x[n:(n - rnum + 1)]
      }
      x = c(x.lmir, x, x.rmir)
      if(reflect.data){
        x = c(x, x[n.ext:1])
      }

      return(list(x = x, idx = (lnum + 1):(lnum + n)))
    }
  }
}


#'Obtain index(position of coefficients in a vector) of each scale.
#'@param t current wavelet level, take value 0:(nt-1), 0 is the coasest level
#'@param nt total number of level
#'@param type 'wavelet', 'station'
get_wave_index = function(t,nt,type){
  # find level index before the current level t
  levels = 0:(nt-1)
  if(t<nt){

    if(type=='wavelet'){
      if(t==levels[nt]){
        return(1:2^t)
      }else{
        idx = levels[levels>t]
        return((sum(2^idx)+1):(sum(2^idx)+2^t))
      }
    }
    if(type=='station'){
      p = 2^nt
      temp = p*(nt-t-1)
      (temp+1):(temp+p)
    }

  }else{
    stop('wavelet levels are from 0 to (#levels-1)')
  }

}


