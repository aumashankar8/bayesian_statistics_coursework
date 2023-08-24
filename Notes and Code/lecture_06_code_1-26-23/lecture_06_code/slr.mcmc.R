slr.mcmc <- function(y,x,n.mcmc){

#
#  Simple Linear Regression
#

####
####  Set up variables
####

n=length(y)

beta.0.save=rep(0,n.mcmc)
beta.1.save=rep(0,n.mcmc)
s2.save=rep(0,n.mcmc)

####
####  Priors (try not to use defaults here!)
####

mu.0=0
s2.0=100
mu.1=0
s2.1=100
q=0.001
r=1000

####
####  Starting Values 
####

beta.0=mean(y)
beta.1=0

####
####  Begin MCMC Loop 
####

for(k in 1:n.mcmc){
  if(k%%1000==0){cat(k," ")}

  ####
  ####  sample s2 
  ####
 
  q.tmp=n/2+q
  r.tmp=1/(sum((y-beta.0-beta.1*x)^2)/2+1/r)
  s2=1/rgamma(1,q.tmp,,r.tmp)

  ####
  ####  sample beta0 
  ####

  var.tmp=1/(n/s2+1/s2.0)
  mn.tmp=var.tmp*(sum(y-beta.1*x)/s2+mu.0/s2.0)
  beta.0=rnorm(1,mn.tmp,sqrt(var.tmp))

  ####
  ####  sample beta1
  ####

  var.tmp=1/(sum(x^2)/s2+1/s2.1)
  mn.tmp=var.tmp*(sum((y-beta.0)*x)/s2+mu.0/s2.0)
  beta.1=rnorm(1,mn.tmp,sqrt(var.tmp))

  ####
  ####  save samples 
  ####

  s2.save[k]=s2
  beta.0.save[k]=beta.0
  beta.1.save[k]=beta.1

};cat("\n")

####
####  Write Output 
####

list(beta.0.save=beta.0.save,beta.1.save=beta.1.save,s2.save=s2.save,n.mcmc=n.mcmc,mu.0=mu.0,s2.0=s2.0,mu.1=mu.1,s2.1=s2.1,q=q,r=r)

}
