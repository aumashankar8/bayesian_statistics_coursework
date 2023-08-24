hw3.betareg.mcmc<-function(y,X,beta.mn,beta.var,beta.tune,tau.gamma1,tau.gamma2,n.mcmc){

#Subroutines
  
logit <- function(theta){
  log(theta/(1-theta))
}
  
logit.inv <- function(z){
  exp(z)/(1+exp(z))
}

#Preliminary variables
n.burn=round(n.mcmc/10)
X=as.matrix(X)
y=as.vector(y)
n=length(y)
p=dim(X)[2]

beta.save=matrix(0,p,n.mcmc)
tau.save=rep(0,n.mcmc)

Davg.save=rep(0,n.mcmc)
mse.y=rep(0,n.mcmc)
mse.ypred=rep(0,n.mcmc)
msediff.save=rep(0,n.mcmc)
mse.save=rep(0,n.mcmc)
ypred.save=matrix(0,n,n.mcmc)
mui.save=matrix(0,n,n.mcmc)

#Starting values
beta=rnorm(p,beta.mn,sqrt(beta.var))
tau=1 #positive number
tau.tune=0.1

mui=logit.inv(X%*%beta)
ai=mui*tau
bi=(1-mui)*tau

beta.acc=1
tau.acc=1

#MCMC Loop
for (k in 1:n.mcmc){
  if(k%%1000==0) cat(k," ")
  
  #Sample Beta
  beta.star=rnorm(p,beta,beta.tune)
  mui.star=logit.inv(X%*%beta.star)
  ai.star=mui.star*tau
  bi.star=(1-mui.star)*tau
  
  mh1=sum(dbeta(y,shape1= ai.star,shape2 = bi.star,log=TRUE))+sum(dnorm(beta.star,beta.mn,sqrt(beta.var),log=TRUE))
  mh2=sum(dbeta(y,shape1= ai,shape2= bi,log=TRUE))+sum(dnorm(beta,beta.mn,sqrt(beta.var),log=TRUE))
  mh.beta=exp(mh1-mh2)
  
  #Updates beta
  if(mh.beta>runif(1)){
    beta=beta.star
    mui=mui.star
    ai=ai.star
    bi=bi.star
    beta.acc=beta.acc+1
  }
  
  #Sample tau
  tau.star=rnorm(1,tau,tau.tune)
  
  #Updates tau
  if(tau.star>0){
    ai.star=mui*tau.star
    bi.star=(1-mui)*tau.star
    
    mh1=sum(dbeta(y,shape1= ai.star,shape2= bi.star,log=TRUE))+sum(dgamma(tau.star,tau.gamma1,tau.gamma2,log=TRUE))
    mh2=sum(dbeta(y,shape1= ai,shape2 = bi,log=TRUE))+sum(dgamma(tau,tau.gamma1,tau.gamma2,log=TRUE))
    mh.tau=exp(mh1-mh2)
    
    if(mh.tau>runif(1)){
      tau=tau.star
      mui=mui.star
      ai=ai.star
      bi=bi.star
      tau.acc=tau.acc+1
    }
    
  }
  
  #DIC Calculations
  
  Davg.save[k]=-2*(sum(dbeta(y,ai,bi,log=TRUE)))

  #Obtain prediction variables
  y.pred=rbeta(n,ai,bi)

  mse.y[k]=mean((y-mui)^2)
  mse.ypred[k]=mean((y.pred-mui)^2)

  msediff.save[k]=mse.ypred[k]-mse.y[k]

  mse.save[k]=mean((y.pred-y)^2)
  # Save Samples
  ypred.save[,k]=y.pred
  beta.save[,k]=beta
  tau.save[k]=tau
}

#Calculate DIC

postmui.mn=apply(mui.save[,-(1:n.burn)],1,mean)
posttau.mn=mean(tau.save[-(1:n.burn)])
D.hat=-2*(sum(dbeta(y,shape1=postmui.mn*posttau.mn,shape2=(1-postmui.mn)*posttau.mn,log=TRUE)))
D.avg=mean(Davg.save[-(1:n.burn)])
pD=D.avg-D.hat
DIC=2*D.avg-D.hat


cat("Dhat:",D.hat,"Davg:",D.avg,"pD:",pD,"DIC:",DIC,"\n")

#Calculate P-value based on MSE
p.value=sum(mse.ypred>mse.y)/n.mcmc

#  Write output
list(y=y,X=X,n.mcmc=n.mcmc,beta.save=beta.save,tau.save=tau.save,beta.acc=beta.acc,p.value=p.value, ypred.save=ypred.save)
}
