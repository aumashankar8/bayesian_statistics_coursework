pois.reg.mcmc <- function(y,X,beta.mn,beta.var,beta.tune,n.mcmc){

###
###  Preliminary Variables
###

n.burn=round(n.mcmc/10)
X=as.matrix(X)
y=as.vector(y)
n=length(y)
p=dim(X)[2]

mse.y=rep(0,n.mcmc)
mse.ypred=rep(0,n.mcmc)
msediff.save=rep(0,n.mcmc)
ypred.save=matrix(0,n,n.mcmc)
mse.save=rep(0,n.mcmc)
beta.save=matrix(0,p,n.mcmc)
lam.save=matrix(0,n,n.mcmc)
Davg=0

###
###  Starting Values
###

beta=c(coef(glm(y~0+X,family=poisson(link="log"))))
beta.sd=sqrt(beta.var)
lam=exp(X%*%beta)
Sig.beta.inv=solve(beta.var*diag(p))

beta.save[,1]=beta
lam.save[,1]=lam
beta.acc=1

###
###  MCMC Loop
###

for(k in 1:n.mcmc){
  if(k%%1000==0) cat(k," ")

  ###
  ### Sample beta
  ###

  beta.star=rnorm(p,beta,beta.tune)
  lam.star=exp(X%*%beta.star)

  mh1=sum(dpois(y,lam.star,log=TRUE))-0.5*t(beta.star-beta.mn)%*%Sig.beta.inv%*%(beta.star-beta.mn)
  mh2=sum(dpois(y,lam,log=TRUE))-0.5*t(beta-beta.mn)%*%Sig.beta.inv%*%(beta-beta.mn)
  mh=exp(mh1-mh2)

  if(mh > runif(1)){
    beta=beta.star   
    lam=lam.star   
    beta.acc=beta.acc+1
  }

  ###
  ### DIC Calculations 
  ###

  if(k>n.burn){
    Davg=Davg-2*(sum(dpois(y,lam,log=TRUE)))/(n.mcmc-n.burn)
  }

  ###
  ### Obtain Predictions 
  ###

  ypred=rpois(n,lam)

  mse.y[k]=mean((y-lam)^2)
  mse.ypred[k]=mean((ypred-lam)^2)
  msediff.save[k]=mse.ypred[k]-mse.y[k]
  mse.save[k]=mean((ypred-y)^2)

  ###
  ### Save Samples 
  ###

  ypred.save[,k]=ypred
  beta.save[,k]=beta
  lam.save[,k]=lam

}
cat("\n")

###
###  Calculate DIC 
###

postlam.mn=apply(lam.save[,-(1:n.burn)],1,mean)
Dhat=-2*(sum(dpois(y,postlam.mn,log=TRUE)))
pD=Davg-Dhat
DIC=2*Davg-Dhat

#cat("Dhat:",Dhat,"Davg:",Davg,"pD:",pD,"DIC:",DIC,"\n")

###
### Calculate P-value based on MSE 
###

p.value=sum(mse.ypred>mse.y)/n.mcmc
 
###
###  Write output 
###

list(y=y,X=X,n.mcmc=n.mcmc,beta.save=beta.save,lam.save=lam.save,beta.acc=beta.acc,msediff.save=msediff.save,mse.y=mse.y,mse.ypred=mse.ypred,ypred.save=ypred.save,mse.save=mse.save,p.value=p.value,DIC=DIC)

}
