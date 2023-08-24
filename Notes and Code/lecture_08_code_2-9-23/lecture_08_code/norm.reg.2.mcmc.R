norm.reg.2.mcmc <- function(y,X,X.pred,beta.mn,beta.var,n.mcmc,no.print=FALSE){

#
#  Multiple Linear Bayesian Regression  (iid errors) w/ DIC calculations.
#

###
### Subroutines 
###

library(mvtnorm)

###
### Hyperpriors
###

n=dim(X)[1]
p=dim(X)[2]
n.pred=dim(X.pred)[1]
n.burn=round(.1*n.mcmc)

r=1000
q=0.001
Sig.beta=beta.var*diag(p)

beta.save=matrix(0,p,n.mcmc)
s2.save=rep(0,n.mcmc)
Dbar.save=rep(0,n.mcmc)
y.pred.save=matrix(0,n.pred,n.mcmc)
y.fit.save=matrix(0,n.pred,n.mcmc)

###
### Starting Values
###

beta=solve(t(X)%*%X)%*%t(X)%*%y

###
### MCMC Loop
###

for(k in 1:n.mcmc){
  if(k%%1000==0) cat(k," ")

  ###
  ### Sample s2
  ###

  tmp.r=(1/r+.5*t(y-X%*%beta)%*%(y-X%*%beta))^(-1)
  tmp.q=n/2+q

  s2=1/rgamma(1,tmp.q,,tmp.r) 

  ###
  ### Sample beta
  ###

  tmp.var=solve(t(X)%*%X/s2 + solve(Sig.beta))
  tmp.mn=tmp.var%*%(t(X)%*%y/s2 + solve(Sig.beta)%*%beta.mn)

  beta=as.vector(rmvnorm(1,tmp.mn,tmp.var,method="chol"))

  ###
  ### DIC Calculations 
  ###

  Dbar.save[k]=-2*sum(dnorm(y,X%*%beta,sqrt(s2),log=TRUE))

  ###
  ### Posterior Predictive Calculations 
  ###

    y.fit.save[,k]=X.pred%*%beta
    y.pred.save[,k]=rnorm(n.pred,y.fit.save[,k],sqrt(s2))

  ###
  ### Save Samples
  ###
  
  beta.save[,k]=beta
  s2.save[k]=s2

}
cat("\n")

###
###  Calculate DIC and Print to Screen
###

if(dim(X)[2]==1){
  postbetamn=mean(beta.save[,-(1:n.burn)])
}
if(dim(X)[2]>1){
  postbetamn=apply(beta.save[,-(1:n.burn)],1,mean)
}
posts2mn=mean(s2.save[-(1:n.burn)])
Dhat=-2*(sum(dnorm(y,X%*%postbetamn,sqrt(posts2mn),log=TRUE)))
Dbar=mean(Dbar.save[-(1:n.burn)])
pD=Dbar-Dhat
DIC=Dhat+2*pD


if(!no.print){
  cat("Posterior Mean for Beta:","\n")
  print(postbetamn)
  cat("Posterior Mean for s2:","\n")
  print(posts2mn)
  cat("Dhat:",Dhat,"Dbar:",Dbar,"pD:",pD,"DIC:",DIC,"\n")
}

###
###  Write Output
###

list(beta.save=beta.save,s2.save=s2.save,y=y,X=X,n.mcmc=n.mcmc,n=n,r=r,q=q,p=p,Dhat=Dhat,Dbar=Dbar,pD=pD,DIC=DIC,y.pred.save=y.pred.save,y.fit.save=y.fit.save)

}
