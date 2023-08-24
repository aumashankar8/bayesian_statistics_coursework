nb.reg.mcmc <- function(y,X,beta.mn,beta.var,beta.tune,logN.mn,logN.sd,n.mcmc){

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
N.save=rep(0,n.mcmc)
beta.save=matrix(0,p,n.mcmc)
lam.save=matrix(0,n,n.mcmc)
Davg.save=rep(0,n.mcmc)

###
###  Starting Values and Priors
###

#logN.mn=log(1)
#logN.sd=log(3)

N=4

#beta=beta.mn  # use this to see good example of burn-in!

beta=c(coef(glm(y~0+X,family=poisson(link="log"))))
lam=exp(X%*%beta)
beta.save[,1]=beta
lam.save[,1]=lam
beta.acc=1
N.tune=.5

###
###  Gibbs Loop
###

for(k in 1:n.mcmc){
  if(k%%1000==0) cat(k," ")

  ###
  ### Sample beta
  ###

  beta.star=rnorm(p,beta,beta.tune)
  lam.star=exp(X%*%beta.star)

  mh1=sum(dnbinom(y,mu=lam.star,size=N,log=TRUE))+sum(dnorm(beta.star,beta.mn,sqrt(beta.var),log=TRUE))
  mh2=sum(dnbinom(y,mu=lam,size=N,log=TRUE))+sum(dnorm(beta,beta.mn,sqrt(beta.var),log=TRUE))

  mh.ratio=exp(mh1-mh2)
  if(mh.ratio > runif(1)){
    beta=beta.star   
    lam=lam.star   
    beta.acc=beta.acc+1
  }

  ###
  ### Sample N 
  ###

  N.star=exp(rnorm(1,log(N),N.tune))

  mh1=sum(dnbinom(y,mu=lam,size=N.star,log=TRUE))+dnorm(log(N.star),logN.mn,logN.sd,log=TRUE)
  mh2=sum(dnbinom(y,mu=lam,size=N,log=TRUE))+dnorm(log(N),logN.mn,logN.sd,log=TRUE)
  
  mhratio=exp(mh1-mh2)
  if(mhratio > runif(1)){
    N=N.star
  }

  ###
  ### DIC Calculations 
  ###

  Davg.save[k]=-2*(sum(dnbinom(y,mu=lam,size=N,log=TRUE)))

  ###
  ### Obtain Predictions 
  ###

  y.pred=rnbinom(n,mu=lam,size=N)

  mse.y[k]=mean((y-lam)^2)
  mse.ypred[k]=mean((y.pred-lam)^2)

  msediff.save[k]=mse.ypred[k]-mse.y[k]
 
  mse.save[k]=mean((y.pred-y)^2)

  ###
  ### Save Samples 
  ###

  ypred.save[,k]=y.pred
  beta.save[,k]=beta
  N.save[k]=N
  lam.save[,k]=lam

}
cat("\n")

###
###  Calculate DIC 
###

postlam.mn=apply(lam.save[,-(1:n.burn)],1,mean)
N.mn=mean(N.save[-(1:n.burn)])
D.hat=-2*(sum(dnbinom(y,mu=postlam.mn,size=N.mn,log=TRUE)))
D.avg=mean(Davg.save[-(1:n.burn)])
pD=D.avg-D.hat
DIC=2*D.avg-D.hat

cat("Dhat:",D.hat,"Davg:",D.avg,"pD:",pD,"DIC:",DIC,"\n")

###
### Calculate P-value based on MSE 
###

p.value=sum(mse.ypred>mse.y)/n.mcmc
 
###
###  Write output 
###

list(y=y,X=X,n.mcmc=n.mcmc,beta.save=beta.save,lam.save=lam.save,beta.acc=beta.acc,msediff.save=msediff.save,mse.y=mse.y,mse.ypred=mse.ypred,ypred.save=ypred.save,mse.save=mse.save,p.value=p.value,N.save=N.save)

}
