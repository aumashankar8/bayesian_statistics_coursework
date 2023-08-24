nb.mcmc <- function(y,mu.mu,s2.mu,mu.N,s2.N,n.mcmc){

###
###  Setup Variables
###

n=length(y)

mu.save=rep(0,n.mcmc)
log.N.save=rep(0,n.mcmc)

###
###  Starting Values 
###

mu=mean(y)
N=1

log.mu=log(mu)
log.N=log(N)

log.mu.tune=.25
log.N.tune=2

###
###  MCMC Loop 
###

for(k in 1:n.mcmc){
  if(k%%1000==0){cat(k," ")}

  ###
  ###  Sample log(N) 
  ###

  log.N.star=rnorm(1,log.N,log.N.tune)
  N.star=exp(log.N.star)
  
  mh.1=sum(dnbinom(y,mu=mu,size=N.star,log=TRUE))+dnorm(log.N.star,mu.N,sqrt(s2.N),log=TRUE)
  mh.2=sum(dnbinom(y,mu=mu,size=N,log=TRUE))+dnorm(log.N,mu.N,sqrt(s2.N),log=TRUE)
  mh=exp(mh.1-mh.2)

  if(mh>runif(1)){
    log.N=log.N.star
    N=N.star
  }

  ###
  ###  Sample mu
  ###

  log.mu.star=rnorm(1,log.mu,log.mu.tune)
  mu.star=exp(log.mu.star)
  
  mh.1=sum(dnbinom(y,mu=mu.star,size=N,log=TRUE))+dnorm(log.mu.star,mu.mu,sqrt(s2.mu),log=TRUE)
  mh.2=sum(dnbinom(y,mu=mu,size=N,log=TRUE))+dnorm(log.mu,mu.mu,sqrt(s2.mu),log=TRUE)
  mh=exp(mh.1-mh.2)

  if(mh>runif(1)){
    log.mu=log.mu.star
    mu=mu.star
  }

  ###
  ###  Save Samples 
  ###

  mu.save[k]=mu
  log.N.save[k]=log.N

};cat("\n")

###
###  Write Output 
###

list(mu.save=mu.save,log.N.save=log.N.save,n.mcmc=n.mcmc)

}
