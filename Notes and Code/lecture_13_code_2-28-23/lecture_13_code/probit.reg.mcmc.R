probit.reg.mcmc <- function(y,X,beta.mn,beta.var,n.mcmc){

###
###  Subroutines
###

rtn <- function(n,mu,sig2,low,high){
  flow=pnorm(low,mu,sqrt(sig2)) 
  fhigh=pnorm(high,mu,sqrt(sig2)) 
  u=runif(n) 
  tmp=flow+u*(fhigh-flow)
  x=qnorm(tmp,mu,sqrt(sig2))
  x
}

###
###  Preliminary Variables
###

n.burn=round(.1*n.mcmc)
X=as.matrix(X)
y=as.vector(y)
n=length(y)
p=dim(X)[2]

beta.save=matrix(0,p,n.mcmc)
z.save=matrix(0,n,n.mcmc)
Davg.save=rep(0,n.mcmc)

y1=(y==1)
y0=(y==0)

###
###  Hyperparameters and Starting Values
###

Sig.beta.inv=solve(beta.var*diag(p))
beta=beta.mn

###
###  Gibbs Loop
###

for(k in 2:n.mcmc){
  if(k%%1000==0) cat(k," ")

  ###
  ### Sample z 
  ###
 
  z=matrix(0,n,1)
  z1=rtn(sum(y),(X%*%beta)[y1],1,0,Inf)
  z0=rtn(sum(1-y),(X%*%beta)[y0],1,-Inf,0)
  z[y1,]=z1
  z[y0,]=z0

  ###
  ### Sample beta
  ###

  tmp.chol=chol(t(X)%*%X + Sig.beta.inv)
  beta=backsolve(tmp.chol,backsolve(tmp.chol,t(X)%*%z + Sig.beta.inv%*%beta.mn,transpose=TRUE)+rnorm(p))

  ###
  ### Deviance calculations 
  ###

  Davg.save[k]=-2*(sum(dbinom(y,1,pnorm(X%*%beta),log=TRUE)))

  ###
  ### Save Samples 
  ###

  beta.save[,k]=beta
  z.save[,k]=z

}
cat("\n")

###
###  Compute DIC 
###

postbetamn=apply(beta.save[,-(1:n.burn)],1,mean)
print("Posterior Mean for Beta:")
print(postbetamn)
D.hat=-2*(sum(dbinom(y,1,pnorm(X%*%postbetamn),log=TRUE)))
D.avg=mean(Davg.save[-(1:n.burn)])
pD=D.avg-D.hat
DIC=2*D.avg-D.hat

cat("Dhat:",D.hat,"Davg:",D.avg,"pD:",pD,"DIC:",DIC,"\n")
 
###
###  Write output 
###

list(y=y,X=X,n.mcmc=n.mcmc,z.save=z.save,beta.save=beta.save,pD=pD,DIC=DIC,D.avg=D.avg,D.hat=D.hat)

}
