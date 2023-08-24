###
###  Subroutines 
###

logit <- function(p){
  log(p)-log(1-p)
}

logit.inv <- function(x){
  exp(x)/(1+exp(x))
}

###
###  Simulate Data
###

n=10
p.true=.7
set.seed(2022)
y=rbinom(n,1,p.true)
y

###
###  Fit Bernoulli Model
###

source("binom.norm.mcmc.R")
n.mcmc=10000
mu=0
s2=2.25

out.1=binom.norm.mcmc(y,mu,s2,1,n.mcmc)
plot(out.1$p.save,type="l")

###
###  Compare with Bernoulli-Beta Model  
###

source("binom.beta.mcmc.R")
n.mcmc=10000
alpha=1
beta=1

out.2=binom.beta.mcmc(y,alpha,beta,n.mcmc)
plot(out.2$p.save,type="l")

plot(density(out.2$p.save),xlim=c(0,1),lwd=2,main="",xlab="p")
lines(density(out.1$p.save),col=2,lwd=2)
curve(dbeta(x,1,1),lty=2,add=TRUE)
lines(density(logit.inv(rnorm(n.mcmc,mu,sqrt(s2)))),col=2,lty=2)
curve(dbeta(x,sum(y)+alpha,sum(1-y)+beta),lwd=2,col=3,add=TRUE)
legend("topleft",lwd=2,col=1:3,legend=c("p MCMC","logit(p) MCMC","p Analytical"))

