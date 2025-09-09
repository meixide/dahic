library(MASS)


generatecov <- function(p,s1, diago=1) {
  A <- matrix(rnorm(p*p,0,s1), p, p)  # Step 1: Generate a random n x n matrix
  A <- t(A) %*% A                # Step 2: Make it symmetric and positive semi-definite
  A <- A + diago * diag(p)       # Step 3: Make the matrix positive definite
  return(A)
}

onenv=function(n,p,betastar,kstar,s1,s2,intercept=0,diago=0.1) {
  Mu=rnorm(p,mean=0,sd=s2)
  Cov=generatecov(p,s1)
  print(Mu)
  d=t(kstar)%*%solve(Cov)%*%kstar + diago  # Schur complement must be positive 
  Covv=rbind(cbind(Cov,kstar),c(kstar,d))
  eps=mvrnorm(n=n,mu=c(Mu,0),Sigma=as.matrix(Covv))
  
  x=eps[,-(p+1)] #+ v
  y=intercept + as.matrix(x)%*%as.vector(betastar) + eps[,p+1]
  return(list(x,y,d))
}

generate_main=function(p,ne,E,s1,s2,betastar,kstar) {




intstar=1

sim=onenv(ne,p,betastar,kstar,s1,s2,intercept=intstar)
x=as.matrix(sim[[1]])
y=sim[[2]]
vary=sim[[3]]
datadrig=list()
Covs=list()
Mus=list()
Covs[[1]]=cov(sim[[1]]) # comment in 1d case
Mus[[1]]=apply(sim[[1]],2,mean) # comment in 1d case

datadrig[[1]]=cbind(y,1,x)

if(E>1) {
  for(i in 2:E) {
    sim=onenv(ne,p,betastar,kstar,s1,s2,intercept=intstar)
    x=rbind(x,as.matrix(sim[[1]]))
    y=rbind(y,sim[[2]])
    vary=rbind(vary,sim[[3]])
    Covs[[i]]=cov(sim[[1]])
    Mus[[i]]=apply(sim[[1]],2,mean)
    datadrig[[i]]=cbind(sim[[2]],1,as.matrix(sim[[1]]))
  }
}


x=as.data.frame(x)
x$ENV=rep(1:E,each=ne)

return(list(x=x,y=y,datadrig=datadrig,Covs=Covs,Mus=Mus,vary=vary,betastar=betastar,kstar=kstar))

}
