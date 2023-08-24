###
###  Setup Simulation 
###

library(mvnfast)

n=2
mu=1
s2=.8
Sig1=s2*diag(n)
Tau1=solve(Sig1)

rho=.9
R=diag(n)
R[1,2]=rho
R[2,1]=rho
Sig2=s2*R
Tau2=solve(Sig2)

###
###  Setup Models
###

library(rjags)
library(coda)

m1.jags <-"
  model{
    y ~ dmnorm(mu.vec,Tau1) 
    mu.vec=rep(mu,n)
    mu ~ dnorm(mu0,tau0)
}
"

m2.jags <-"
  model{
    y ~ dmnorm(mu.vec,Tau2) 
    mu.vec=rep(mu,n)
    mu ~ dnorm(mu0,tau0)
}
"

n.mcmc=10000
n.burn=round(.2*n.mcmc)
mod1<-textConnection(m1.jags)
mod2<-textConnection(m2.jags)

###
###  Simulate Single Dataset and Fit Models using JAGS
###

mu0=0
tau0=1/1000

set.seed(404)
y=c(rmvn(1,rep(mu,2),Sig2)) # simulate data from model w/ dependence

m1.out<-jags.model(mod1,data=list('y'=y,'n'=n,'mu0'=mu0,'tau0'=tau0,'Tau1'=Tau1),n.chains=1,n.adapt=0) 
m2.out<-jags.model(mod2,data=list('y'=y,'n'=n,'mu0'=mu0,'tau0'=tau0,'Tau2'=Tau2),n.chains=1,n.adapt=0) 

m1.samp=jags.samples(m1.out,c('mu'),n.mcmc)
p1=mean(m1.samp$mu[1,,1]>0)
m2.samp=jags.samples(m2.out,c('mu'),n.mcmc) 
p2=mean(m2.samp$mu[1,,1]>0)

p1
p2

d1=density(m1.samp$mu[1,,1])
d2=density(m2.samp$mu[1,,1])

plot(d1,type="l",lwd=2,col=rgb(0,0,0,.5),ylim=c(0,max(d1$y,d2$y)),main="",xlab=bquote(mu))
lines(d2,type="l",lwd=2,col=rgb(1,0,0,.5))
abline(v=0,col=rgb(0,0,0,.4),lty=2)
legend("topright",col=c(rgb(0,0,0,.5),rgb(1,0,0,.5)),lwd=2,legend=c("w/o corr","w/ corr"))

