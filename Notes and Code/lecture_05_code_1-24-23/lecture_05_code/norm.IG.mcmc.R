norm.IG.mcmc <- function(y,mu.0,s2.0,mu.strt,n.mcmc){

#
#  Gibbs Sampler for Gaussian Data Model, with Gaussian prior on mu and IG on s2.
#
#  Example Use:
#  mcmc.out=norm.IG.mcmc(rnorm(10,3,2),0,100,0,10000)
#  mcmc.out=norm.IG.mcmc(rnorm(100,3,2),0,100,0,10000)
#

####
####  Set up variables 
####

n.burn=round(n.mcmc/10)
n=length(y)
mu.save=rep(0,n.mcmc)
s2.save=rep(0,n.mcmc)
mu=mu.strt

r=1000
q=0.001

####
####  Begin Gibbs Loop 
####

for(k in 1:n.mcmc){
  if((k%%1000)==0) cat(k," ")

  ####
  ####  Sample s2 
  ####
  
  tmp.r=1/(sum((y-mu)^2)/2+1/r)
  tmp.q=n/2+q
  
  s2=1/rgamma(1,tmp.q,,tmp.r)

  ####
  ####  Sample mu 
  ####

  tmp.mn=(s2*mu.0+s2.0*sum(y))/(s2+n*s2.0)
  tmp.var=s2*s2.0/(s2+n*s2.0)

  mu=rnorm(1,tmp.mn,sqrt(tmp.var))

  ####
  ####  Save Samples 
  ####

  mu.save[k]=mu
  s2.save[k]=s2

}
cat("\n")

####
####  Write Output 
####

list(mu.save=mu.save,s2.save=s2.save,mu.0=mu.0,s2.0=s2.0,q=q,r=r,n.mcmc=n.mcmc,n.burn=n.burn)

}
