###
###  Simulate Data
###

set.seed(2022)
n=100
lam=5
y=rpois(n,lam)

plot(table(y)/n,type="h",ylab="probability")

###
###  Fit Poisson Model 
###

lam.mn=10
lam.var=10^2

gam.1=(lam.mn^2)/lam.var
gam.2=lam.mn/lam.var

n.mcmc=10000
out.1=vector("list")
out.1$lam.save=rgamma(n.mcmc,sum(y)+gam.1,n+gam.2)

###
###  Fit Negative Binomial Model 
###

source("nb.mcmc.R")
n.mcmc=10000

mu.mu=10
s2.mu=100
mu.N=0
s2.N=10

out.2=nb.mcmc(y,mu.mu,s2.mu,mu.N,s2.N,n.mcmc)

layout(matrix(1:2,2,1))
plot(out.2$mu.save,type="l")
plot(out.2$log.N.save,type="l")

###
###  Compare Posterior Results 
###

layout(matrix(1:2,1,2))
plot(density(out.1$lam.save),col=rgb(0,0,0,.5),lwd=2,main="",xlab=bquote(mu))
lines(density(out.2$mu.save),col=rgb(1,0,0,.5),lwd=2)
plot(density(out.2$log.N.save),col=rgb(1,0,0,.5),lwd=2,main="",xlab="log(N)")

###
###  Compare DIC 
###

D.hat.1=-2*sum(dpois(y,mean(out.1$lam.save),log=TRUE))
D.bar.1=0
for(k in 1:n.mcmc){
  D.bar.1=D.bar.1-2*sum(dpois(y,out.1$lam.save[k],log=TRUE))/n.mcmc
}
p.D.1=D.bar.1-D.hat.1
DIC.1=D.hat.1+2*p.D.1
c(p.D.1,DIC.1)

D.hat.2=-2*sum(dnbinom(y,mu=mean(out.2$mu.save),size=mean(exp(out.2$log.N.save)),log=TRUE))
D.bar.2=0
for(k in 1:n.mcmc){
  D.bar.2=D.bar.2-2*sum(dnbinom(y,mu=out.2$mu.save[k],size=exp(out.2$log.N.save[k]),log=TRUE))/n.mcmc
}
p.D.2=D.bar.2-D.hat.2
DIC.2=D.hat.2+2*p.D.2
c(p.D.2,DIC.2)





