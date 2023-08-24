binom.beta.mcmc <- function(y,alpha,beta,n.mcmc){

####
#######################
####
####  Implements M-H sampler for binomial-beta, where y~Binom(1,p), p~beta(alpha,beta).
####
####  y: vector of binomial count data of length n
####  p: probability of success
####  alpha, beta: hyperparameters in prior 
####
####  Example Use:
####  
####  binom.beta.mcmc(c(1,0,1,1,1,0,1,1,1,0),1,1,1000) 
####  
#######################

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

####
####  Begin Gibbs Loop 
####
  
for(k in 1:n.mcmc){
  if(k%%1000==0)  cat(k," ")

  ####
  ####  Sample phi 
  ####

  p.star=runif(1)

  mh1=sum(dbinom(y,1,p.star,log=TRUE))+dbeta(p.star,alpha,beta,log=TRUE)
  mh2=sum(dbinom(y,1,p,log=TRUE))+dbeta(p,alpha,beta,log=TRUE)
  mh=exp(mh1-mh2)
  
 if(mh > runif(1)){
    p=p.star
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
