hd_weighted=function(Y,X,intercept=T,weights) {
  # never used, does not make sense
  if(intercept) {
    X=cbind(1,X)
  }
  M <- X %>%
    group_by(ENV) %>%
    mutate(across(everything(), \(x) mean(x,na.rm = TRUE))) %>%
    ungroup() %>%
    dplyr::select(-ENV) %>%
    as.matrix()
  Y=weights*Y
  X=as.data.frame(as.matrix(diag(weights))%*%as.matrix(X))
  mminv=solve(t(M)%*%M)
  betahat=mminv%*%t(M)%*%Y
  
  if(intercept) {
    x=X[,-c(1,ncol(X))]
  }else {
    x=X[,-c(ncol(X))]
  }
  
  n=nrow(X)
  
  E=length(table(X$ENV))
  
  if(intercept) {
    m=as.matrix(M[,-1])
    
  }else {
    m=as.matrix(M)
  }
  
  x=as.matrix(x)
  khat=as.double(c(1/(n))*t(x)%*%(diag(1,n) - x%*%solve(t(m)%*%m)%*%t(m))%*%Y)
  fac=solve(t(x)%*%x - t(m)%*%m)
  
  kopt=fac%*%t(x)%*%(diag(1,n) - x%*%solve(t(m)%*%m)%*%t(m))%*%Y
  if(!intercept) {
    condmean=x%*%betahat
  }else{condmean=betahat[1]+x%*%betahat[-1]}
  vary=as.double(mean((Y-condmean-(x-m)%*%kopt
  )^2) + sum(kopt*khat))
  
  return(list(betahat=betahat,khat=khat,kopt=kopt,vary=vary))
  
}

library('dplyr')
hd=function(Y,X,intercept=T) {
  if(intercept) {
    X=cbind(1,X)
  }
  M <- X %>%
    group_by(ENV) %>%
    mutate(across(everything(), \(x) mean(x,na.rm = TRUE))) %>%
    ungroup() %>%
    dplyr::select(-ENV) %>%
    as.matrix()
  
  mminv=solve(t(M)%*%M)
  betahat=mminv%*%t(M)%*%Y
  
  if(intercept) {
    x=X[,-c(1,ncol(X))]
  }else {
    x=X[,-c(ncol(X))]
  }

  n=nrow(X)
  
  E=length(table(X$ENV))
  
  if(intercept) {
    m=as.matrix(M[,-1])
    
  }else {
    m=as.matrix(M)
  }
  
  x=as.matrix(x)
  khat=as.double(c(1/(n-E))*t(x)%*%(diag(1,n) - x%*%solve(t(m)%*%m)%*%t(m))%*%Y)
  fac=solve(t(x)%*%x - t(m)%*%m)
  
  kopt=fac%*%t(x)%*%(diag(1,n) - x%*%solve(t(m)%*%m)%*%t(m))%*%Y
  if(!intercept) {
    condmean=x%*%betahat
  }else{condmean=betahat[1]+x%*%betahat[-1]}
  vary=as.double(mean((Y-condmean-(x-m)%*%kopt
  )^2) + sum(kopt*khat))
  
  return(list(betahat=betahat,khat=khat,kopt=kopt,vary=vary))
  
}

generator=function(betahat,khat,xnew,vary,distributional=F) {
  epsyprime=as.matrix(sweep(xnew, 2, apply(xnew,2,mean), FUN = "-"))%*%solve(cov(xnew))%*%khat
  #print(cov(epsyprime,xnew))
  #print(var(epsyprime))
  noise=0#rnorm(nrow(xnew),0,sqrt(vary-var(epsyprime)))
  #var(noise)
  return(epsyprime+betahat[1]+as.matrix(xnew)%*%betahat[-1]+distributional*noise)
}

