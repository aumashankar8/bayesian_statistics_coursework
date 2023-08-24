###
### Simulate Data
###

set.seed(1001)
n=10
theta.true=.7
y=rbinom(n,1,theta.true)

y

###
###  Fit Model using MCMC 
###

n.mcmc=200000
source("bern.rw.mcmc.R")
out.rw=bern.rw.mcmc(y,.5,n.mcmc)
source("bern.tn.mcmc.R")
out.tn=bern.tn.mcmc(y,.5,n.mcmc)

###
###  Compare Results 
###

plot(density(out.rw$theta.save),type="l",col=1,lwd=2,xlim=c(0,1),main="",xlab=bquote(theta))
lines(density(out.tn$theta.save),type="l",col=2,lwd=2)
curve(dbeta(x,sum(y)+out.rw$alpha,sum(1-y)+out.rw$beta),add=TRUE,col=3,lty=2,lwd=2)

