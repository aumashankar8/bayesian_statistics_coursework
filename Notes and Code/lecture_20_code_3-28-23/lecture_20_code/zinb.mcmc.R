zinb.mcmc <- function(y,n.mcmc){

#
#  ZINB model with homogeneous lambda
#

###
###  Set up variables
###

n=length(y)

lam.save=rep(0,n.mcmc)
p.save=rep(0,n.mcmc)
N.save=rep(0,n.mcmc)
z.mean=rep(0,n)

n.burn=round(.2*n.mcmc)

###
###  Priors 
###

alpha.p=1
beta.p=1

mu.loglam=0
s2.loglam=100

mu.logN=0
s2.logN=1000

###
###  Starting values 
###

lam=mean(y)
loglam=log(lam)
N=1
logN=log(1)
z=rep(0,n)
z[y>0]=1

ll.tune=.1
lN.tune=.1

###
###  Begin MCMC loop 
###

for(k in 1:n.mcmc){
  if(k%%1000==0){cat(k," ")}

  ###
  ###  Sample p
  ###

  p=rbeta(1,sum(z)+alpha.p,sum(1-z)+beta.p)

  ###
  ###  Sample log(lambda)
  ###
  
  loglam.star=rnorm(1,loglam,ll.tune)
  mh.1=sum(dnbinom(y[z==1],mu=exp(loglam.star),size=N,log=TRUE))+dnorm(loglam.star,mu.loglam,sqrt(s2.loglam),log=TRUE)
  mh.2=sum(dnbinom(y[z==1],mu=lam,size=N,log=TRUE))+dnorm(loglam,mu.loglam,sqrt(s2.loglam),log=TRUE)
  mh=exp(mh.1-mh.2)
  if(mh>runif(1)){
    loglam=loglam.star 
    lam=exp(loglam)
  }

  ###
  ###  Sample log(N)
  ###
  
  logN.star=rnorm(1,logN,lN.tune)
  mh.1=sum(dnbinom(y[z==1],mu=lam,size=exp(logN.star),log=TRUE))+dnorm(logN.star,mu.logN,sqrt(s2.logN),log=TRUE)
  mh.2=sum(dnbinom(y[z==1],mu=lam,size=N,log=TRUE))+dnorm(logN,mu.logN,sqrt(s2.logN),log=TRUE)
  mh=exp(mh.1-mh.2)
  if(mh>runif(1)){
    logN=logN.star 
    N=exp(logN)
  }

  ###
  ###  Sample z 
  ###

  p.tmp=p*(N/(N+lam))^N/(p*(N/(N+lam))^N+1-p)
  z[y==0]=rbinom(sum(y==0),1,p.tmp)

  ###
  ###  Save Samples
  ###

  p.save[k]=p
  N.save[k]=N
  lam.save[k]=lam
  if(k>n.burn){
    z.mean=z.mean+z/(n.mcmc-n.burn)
  }

};cat("\n")

###
### Write Output 
###

list(p.save=p.save,lam.save=lam.save,N.save=N.save,z.mean=z.mean,n.mcmc=n.mcmc)

}
