mix.mcmc <- function(y,n.mcmc){

#
#  Fits mixture model with 2 Gaussian components
#

####
####  Setup Variables
####

n=length(y)
mu.save=matrix(0,n.mcmc,2)
s2.save=matrix(0,n.mcmc,2)
p.save=rep(0,n.mcmc)
z.mean=rep(0,n)
n.burn=round(.2*n.mcmc)

####
####  Priors 
####

mu.10=70
s2.10=10000
mu.20=70
s2.20=10000

q.1=.01
r.1=100
q.2=.01
r.2=100

alpha=1
beta=1

####
####  Starting Values 
####

p=.5
z=rep(0,n)
z[y<mean(y)]=1
mu.1=mean(y)-20
mu.2=mean(y)+20
s2.1=var(y)/2
s.1=sqrt(s2.1)
s2.2=var(y)/2
s.2=sqrt(s2.2)

####
####  Begin MCMC Loop
####

for(k in 1:n.mcmc){
  if(k%%1000==0) cat(k," ")

  ####
  ####  Sample z
  ####

  p.tmp=p*dnorm(y,mu.1,s.1)/(p*dnorm(y,mu.1,s.1)+(1-p)*dnorm(y,mu.2,s.2))
  z=rbinom(n,1,p.tmp)

  ####
  ####  Sample mu.1 
  ####

  s2.tmp=(sum(z)/s2.1 + 1/s2.10)^(-1)
  mu.tmp=s2.tmp*(sum(y[z==1])/s2.1 + mu.10/s2.10)
  mu.1=rnorm(1,mu.tmp,sqrt(s2.tmp))

  ####
  ####  Sample mu.2
  ####

  s2.tmp=(sum(1-z)/s2.2 + 1/s2.20)^(-1)
  mu.tmp=s2.tmp*(sum(y[z==0])/s2.2+mu.20/s2.20)
  mu.2=rnorm(1,mu.tmp,sqrt(s2.tmp))

  ####
  ####  Sample s2.1 
  ####
  
  q.tmp=sum(z)+q.1
  r.tmp=1/(sum((y[z==1]-mu.1)^2)/2+1/r.1)
  s2.1=1/rgamma(1,q.tmp,,r.tmp)
  s.1=sqrt(s2.1)

  ####
  ####  Sample s2.2 
  ####
  
  q.tmp=sum(1-z)+q.2
  r.tmp=1/(sum((y[z==0]-mu.2)^2)/2+1/r.2)
  s2.2=1/rgamma(1,q.tmp,,r.tmp)
  s.2=sqrt(s2.2)

  ####
  ####  Sample p 
  ####
 
  p=rbeta(1,sum(z)+alpha,sum(1-z)+beta)
 
  ####
  ####  Save Samples 
  ####

  mu.save[k,]=c(mu.1,mu.2)
  s2.save[k,]=c(s2.1,s2.2)
  p.save[k]=p

  if(k>n.burn){
    z.mean=z.mean+z/(n.mcmc-n.burn) 
  }

}
cat("\n")

####
####  Write Output 
####

list(mu.save=mu.save,p.save=p.save,s2.save=s2.save,z.mean=z.mean,n.mcmc=n.mcmc,n.burn=n.burn)

}
