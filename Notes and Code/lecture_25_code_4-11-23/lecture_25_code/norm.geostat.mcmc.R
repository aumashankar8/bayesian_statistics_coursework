norm.geostat.mcmc <- function(y,X,locs,beta.mn,beta.var,n.mcmc){

###
### Subroutines 
###

library(fields)

ldmvnorm<-function(y,mean,Sig.inv,logdet){ 
    logretval=-(logdet+t(y-mean)%*%Sig.inv%*%(y-mean))/2
    logretval
}

###
### Set up Variables and Hyperparameters
###

n=dim(X)[1]
p=dim(X)[2]
r=1000
q=0.001
Sig.beta=beta.var*diag(p)
Sig.beta.inv=solve(Sig.beta)

D=rdist.earth(locs,locs)
max.d=max(D)
phi.mn=max.d/4/4
phi.var=10*phi.mn
gamma.1=(phi.mn^2)/phi.var
gamma.2=phi.mn/phi.var
phi=phi.mn
phi.tune=phi/4

beta.save=matrix(0,p,n.mcmc)
s2.save=rep(0,n.mcmc)
phi.save=rep(0,n.mcmc)

###
### Starting Values
###

beta=solve(t(X)%*%X)%*%t(X)%*%y
Xbeta=X%*%beta
s2=mean((y-Xbeta)^2)

R=exp(-D/phi)
R.inv=solve(R)
logdet=determinant(R)$modulus

###
### MCMC Loop
###

for(k in 1:n.mcmc){
  if(k%%1000==0) cat(k," ")

  ###
  ### Sample s2
  ###

  tmp.r=(1/r+.5*t(y-Xbeta)%*%R.inv%*%(y-Xbeta))^(-1)
  tmp.q=n/2+q
  s2=1/rgamma(1,tmp.q,,tmp.r) 
  Sig.inv=R.inv/s2

  ###
  ### Sample beta
  ###

  tmp.chol=chol(t(X)%*%Sig.inv%*%X + Sig.beta.inv)
  beta=backsolve(tmp.chol,backsolve(tmp.chol,t(X)%*%Sig.inv%*%y + Sig.beta.inv%*%beta.mn,transpose=TRUE)+rnorm(p))
  Xbeta=X%*%beta

  ###
  ### Sample phi 
  ###

  phi.star=rnorm(1,phi,phi.tune)
  if(phi.star>0){ 
    R.star=exp(-D/phi.star)
    R.inv.star=solve(R.star)
    logdet.star=determinant(R.star)$modulus
    mh.1=ldmvnorm(y,Xbeta,R.inv.star/s2,logdet.star)+dgamma(phi.star,gamma.1,gamma.2,log=TRUE)
    mh.2=ldmvnorm(y,Xbeta,R.inv/s2,logdet)+dgamma(phi,gamma.1,gamma.2,log=TRUE)
    if(exp(mh.1-mh.2)>runif(1)){
      phi=phi.star
      R.inv=R.inv.star
      logdet=logdet.star
    }
  }

  ###
  ### Save Samples
  ###
  
  beta.save[,k]=beta
  s2.save[k]=s2
  phi.save[k]=phi

}
cat("\n")

###
###  Write Output
###

list(beta.save=beta.save,s2.save=s2.save,phi.save=phi.save,y=y,X=X,n.mcmc=n.mcmc,n=n,r=r,q=q,p=p,gamma.1=gamma.1,gamma.2=gamma.2,phi.tune=phi.tune)

}
