hier.mcmc <- function(y,x,ind,n.mcmc){

#
#  Fits hierarchical model with random coefficients and single covariate
#

###
###  Load Packages 
###

library(mvnfast)

###
###  Setup Variables
###

n=max(ind)
J=sum(ind==1) # assumes all J homog for all ind

X.list=vector("list",n)

beta.save=array(0,c(2,n,n.mcmc))
mu.save=matrix(0,2,n.mcmc)
rho.save=rep(0,n.mcmc)
s2.beta.save=rep(0,n.mcmc)
s2.save=rep(0,n.mcmc)

###
###  Priors 
###

q=.001
r=1000

mu.0=rep(0,2)
s2.0=1000
Sig.0=s2.0*diag(2)
Sig.0.inv=solve(Sig.0)

q.beta=.001
r.beta=1000

rho.l=-1
rho.u=1

###
###  Starting Values 
###

beta=matrix(0,2,n)

for(i in 1:n){
  X.list[[i]]=matrix(1,J,2)  
  X.list[[i]][,2]=x[ind==i]
  X.tmp=X.list[[i]]
  y.tmp=y[ind==i]
  beta[,i]=solve(t(X.tmp)%*%X.tmp)%*%t(X.tmp)%*%y.tmp
}

mu=apply(beta,1,mean)
Sig.emp=cov(t(beta-mu))
s2.beta=mean(diag(Sig.emp))
rho=cor(t(beta-mu))[1,2]
Sig=diag(2)
Sig[1,2]=rho
Sig[2,1]=rho
Sig=s2.beta*Sig
Sig.inv=solve(Sig)

###
###  Begin MCMC Loop 
###

for(k in 1:n.mcmc){
  if(k%%1000==0){cat(k," ")}

  ###
  ###  Sample s2
  ###

  tmp.sum=0
  for(i in 1:n){
    tmp.sum=tmp.sum+sum((y[ind==i]-X.list[[i]]%*%beta[,i])^2)
  }
  r.tmp=1/(tmp.sum/2+1/r)
  q.tmp=n*J+q
  s2=1/rgamma(1,q.tmp,,r.tmp)

  ###
  ###  Sample beta 
  ###

  for(i in 1:n){
    tmp.cov=solve(t(X.list[[i]])%*%X.list[[i]]/s2+Sig.inv)
    tmp.mn=tmp.cov%*%(t(X.list[[i]])%*%y[ind==i]/s2+Sig.inv%*%mu)
    beta[,i]=t(rmvn(1,tmp.mn,tmp.cov))
  }   

  ###
  ###  Sample mu 
  ###

  tmp.cov=solve(n*Sig.inv+Sig.0.inv)
  tmp.mn=tmp.cov%*%(Sig.inv%*%apply(beta,1,sum)+Sig.0.inv%*%mu.0)
  mu=t(rmvn(1,tmp.mn,tmp.cov))

  ###
  ###  Sample s2.beta 
  ###

  R.inv=Sig.inv*s2.beta
  tmp.sum=0
  for(i in 1:n){
    tmp.sum=tmp.sum+(t(beta[,i]-mu)%*%R.inv%*%(beta[,i]-mu))
  }
  r.tmp=1/(tmp.sum/2+1/r.beta)
  q.tmp=n+q.beta
  s2.beta=1/rgamma(1,q.tmp,,r.tmp)
  Sig.inv=R.inv/s2.beta

  ###
  ###  Sample rho 
  ###

  rho.star=runif(1,rho.l,rho.u) 
  R.star=diag(2)
  R.star[1,2]=rho.star
  R.star[2,1]=rho.star
  Sig.star=s2.beta*R.star
 
  tmp.sum=0
  tmp.sum.star=0
 
  for(i in 1:n){
    tmp.sum=tmp.sum+dmvn(beta[,i],mu,Sig,log=TRUE)
    tmp.sum.star=tmp.sum.star+dmvn(beta[,i],mu,Sig.star,log=TRUE)
  }

  mh.1=tmp.sum.star
  mh.2=tmp.sum
  mh=exp(mh.1-mh.2)
  if(mh>runif(1)){
    rho=rho.star
    Sig=Sig.star
    Sig.inv=solve(Sig)
  }

  ###
  ###  Save Samples 
  ###

  beta.save[,,k]=beta
  mu.save[,k]=mu
  s2.save[k]=s2
  s2.beta.save[k]=s2.beta
  rho.save[k]=rho

};cat("\n")

###
###  Write Output 
###

list(beta.save=beta.save,mu.save=mu.save,s2.save=s2.save,s2.beta.save=s2.beta.save,rho.save=rho.save,n.mcmc=n.mcmc)

}
