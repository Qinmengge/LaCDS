
LaCDS=function(t,ind,X,G=2,tol=0.0001,max_iter = 1000, min_iter=1,alpha=0.5){
  rst=list()
  AIC=c()

  maxt=max(t)
  c=ncol(X)
  r=nrow(X)

  # order t
  od=order(t,decreasing = T)
  t=t[od]
  ind=ind[od]
  X=as.matrix(X[od,])

  #kernal
  nker=maxt
  rho=100
  kernal=rep(seq(1,maxt,length.out=nker),each=maxt)
  #kernal=rep(c(9,17,23,29,34,39,44,49,52,53),each=maxt)
  dim(kernal)=c(maxt,nker)
  t_seq=rep(1:maxt,nker)
  dim(t_seq)=c(maxt,nker)
  wker=exp(-rho*(t_seq-kernal)^2)
  wker=wker/rowSums(wker)

  #formatting
  ind=t(t(ind))
  ker=t(t(c(rep(0,nker))))
  beta_t=t(t(c(rep(0,maxt))))
  beta_v=t(t(rep(0,c)))
  t=t(t(t))

  for (g in 1:G) {
    z=sample(1:g,r,replace = T)
    Z=matrix(rep(0,r*g),nrow = r)
    for (i in 1:g){
      Z[z==i,i]=1
    }
    res=EM(wker,Z=Z,group=g ,t, X, ind, beta_t, beta_v, tol = 0.0001, tol2=0.0001, max_iter = 1000, min_iter=0, alpha=0.5)
    lik1=lik(t, X, ind, res$Beta_t, res$Beta_v, t(t(colMeans(res$Z))), maxt, r, c, g)
    AIC=(length(res$Beta_t)+length(res$Beta_v))*2-2*lik1$likeli2
    sd=sqrt(diag(solve(res$Info+diag(rep(0.0000001,maxt+c)))))
    res$sd=sd
    res$AIC=AIC
    rst[[as.character(g)]]=res
  }
  return(rst)
}
