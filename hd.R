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

library('dplyr')


est_drig <- function(data, gamma, y_idx=1, del_idx=NULL, unif_weight=FALSE) {
  if (is.null(del_idx)) {
    del_idx <- y_idx
  }
  m <- length(data)
  if (unif_weight) {
    w <- rep(1/m, m)
  } else {
    w <- sapply(data, function(e) nrow(e))
    w <- w / sum(w)
  }
  
  gram_x <- list()
  gram_xy <- list()
  
  for (e in 1:m) {
    data_e <- data[[e]]
    n <- nrow(data_e)
    y <- data_e[, y_idx]
    x <- data_e[, -c(y_idx, del_idx)]
    
    gram_x[[e]] <- t(x) %*% x / n
    gram_xy[[e]] <- t(x) %*% y / n
  }
  G <- (1 - gamma) * gram_x[[1]] + gamma * Reduce('+', Map('*', gram_x, w))
  
  
  
  Z <- (1 - gamma) * gram_xy[[1]] + gamma * Reduce('+', Map('*', gram_xy, w))
  drig= solve(G) %*% Z

  if(is.infinite(gamma)) {
    G <-  gram_x[[1]] - gram_x[[2]]
    Z <- gram_xy[[1]] - gram_xy[[2]]
    dantzig=solve(G) %*% Z
    return(dantzig)
  }
  return(drig)
}


