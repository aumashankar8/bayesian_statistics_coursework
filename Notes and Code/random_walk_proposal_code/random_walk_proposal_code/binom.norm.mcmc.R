binom.norm.mcmc <- function(y,mu,s2,sig.tune,n.mcmc){

####
#######################
####
####  Implements M-H sampler for binomial-beta, where y~Binom(1,p), logit(p)~N(mu,s2).
####
####  y: vector of binomial count data of length n
####  p: probability of success
####  mu, s2: hyperparameters in prior 
####  sig.tune:  M-H random walk tuning parameter
####
####  Example Use:
####  
####  binom.norm.mcmc(c(1,0,1,1,1,0,1,1,1,0),0,2.25,1000) 
####  
#######################

####
####  Subroutine 
####

logit <- function(p){
  log(p)-log(1-p)
}

logit.inv <- function(x){
  exp(x)/(1+exp(x)) 
}

####
####  Setup Variables 
####

n=length(y)
p.save=rep(0,n.mcmc)
mh.p=1

####
####  Priors and Starting Values 
####

p=mean(y)
logit.p=logit(p)
sig=sqrt(s2)

####
####  Begin Gibbs Loop 
####
  
for(k in 1:n.mcmc){
  if(k%%1000==0)  cat(k," ")

  ####
  ####  Sample logit(p) 
  ####

  logit.p.star=rnorm(1,logit.p,sig.tune)
  p.star=logit.inv(logit.p.star)

  mh1=sum(dbinom(y,1,p.star,log=TRUE))+dnorm(logit.p.star,mu,sig,log=TRUE)
  mh2=sum(dbinom(y,1,p,log=TRUE))+dnorm(logit.p,mu,sig,log=TRUE)
  mh=exp(mh1-mh2)
  
  if(mh > runif(1)){
    p=p.star
    logit.p=logit.p.star
    mh.p=mh.p+1
  }

  ####
  ####  Save Samples 
  ####

  p.save[k]=p

}
cat("\n")

####
####  Write Output 
####
 
list(p.save=p.save,mh.p=mh.p,n.mcmc=n.mcmc)

}
