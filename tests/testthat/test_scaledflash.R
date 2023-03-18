set.seed(123)
n = 99
p = 300
k= 3
L = matrix(0, nrow=n, ncol=k)
F = matrix(0, nrow=p, ncol=k)

L[1:(n/3),1] = 1
L[((n/3)+1):(2*n/3),2] = 1
L[((2*n/3)+1):n,3] = 1

F[1:(p/3),1] = 1+10*runif(p/3)
F[((p/3)+1):(2*p/3),2] = 1+10*runif(p/3)
F[((2*p/3)+1):p,3] = 1+10*runif(p/3)

lambda = L %*% t(F)
X = matrix(rpois(n=length(lambda),lambda),nrow=n)
image(X)

fit_tm = fastTopics::fit_poisson_nmf(X,3)
plot(fit_tm$L[,1],main="estimated loadings 1")
plot(fit_tm$L[,2],main="estimated loadings 2")
plot(fit_tm$L[,3],main="estimated loadings 3")

Y_tilde = biwhitening(X)
fit_sf = scaledflash(Y_tilde$Y,Y_tilde$u,Y_tilde$v,
                     S2 = NULL,
                     var.type = 'by_column',
                     Kmax=5,
                     tol=0.01,
                     maxiter = 1000,
                     ebnm_fn = 'ebnm_pe',
                     init_fn = 'nnmf_r1',
                     ebnm_param=NULL,
                     verbose=TRUE,
                     nullcheck=TRUE,
                     sigma2 = NULL,
                     seed=12345)
plot(fit_sf$ldf$f[,1])
plot(fit_sf$ldf$f[,2])
plot(fit_sf$ldf$f[,3])




#################
set.seed(123)
n = 99
p = 300
k= 4
mfac = 2 # controls PVE of dense factor
L = matrix(0, nrow=n, ncol=k)
F = matrix(0, nrow=p, ncol=k)

L[1:(n/3),1] = 1
L[((n/3)+1):(2*n/3),2] = 1
L[((2*n/3)+1):n,3] = 1
L[,4] = 1+mfac*runif(n)

F[1:(p/3),1] = 1+10*runif(p/3)
F[((p/3)+1):(2*p/3),2] = 1+10*runif(p/3)
F[((2*p/3)+1):p,3] = 1+10*runif(p/3)
F[,4]= 1+mfac*runif(p)

lambda = L %*% t(F)
X = matrix(rpois(n=length(lambda),lambda),nrow=n)
image(X)

Y_tilde = biwhitening(X)
fit_sf = scaledflash(Y_tilde$Y,Y_tilde$u,Y_tilde$v,
                     S2 = NULL,
                     var.type = 'by_column',
                     Kmax=10,
                     tol=0.01,
                     maxiter = 1000,
                     ebnm_fn = 'ebnm_pe',
                     init_fn = 'nnmf_r1',
                     ebnm_param=NULL,
                     verbose=TRUE,
                     nullcheck=TRUE,
                     sigma2 = NULL,
                     seed=12345)
plot(fit_sf$ldf$l[,1])
plot(fit_sf$ldf$l[,2])
plot(fit_sf$ldf$l[,3])
plot(fit_sf$ldf$l[,4])
