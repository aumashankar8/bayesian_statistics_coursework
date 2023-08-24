###
###  Simulate Data
###

set.seed(1001)

n=100
lam=5
p=.8
z=rbinom(n,1,p)
y=rep(0,n)
y[z==1]=rpois(sum(z),lam)

plot(table(y)/n,type="h")

###
### Fit Model 
###

source("zip.mcmc.R")
n.mcmc=10000

out=zip.mcmc(y,n.mcmc)

layout(matrix(1:2,2,1))
plot(out$lam.save,type="l")
abline(h=lam,col=rgb(1,0,0,.7))
plot(out$p.save,type="l")
abline(h=p,col=rgb(1,0,0,.7))

###
### Posterior Predictive Distribution 
###

z.pred=rbinom(n.mcmc,1,out$p.save)
y.pred=rep(0,n.mcmc)
y.pred[z.pred==1]=rpois(sum(z.pred),out$lam.save)

plot(table(y)/n,type="h")
lines(table(y.pred)/n.mcmc,type="h",lty=2,lwd=4,col=rgb(1,0,0,.5))

###
###  Read in Doctor Visit Data
###

dv.df=read.csv("count_docvisit.csv")
y=dv.df$docvis
plot(table(y),type="h")

###
###  Fit Homogeneous ZIP Model
###

source("zip.mcmc.R")
n.mcmc=10000
out.1=zip.mcmc(y,n.mcmc)

layout(matrix(1:2,2,1))
plot(out.1$lam.save,type="l")
plot(out.1$p.save,type="l")

###
### Posterior Predictive Distribution
###

z.pred=rbinom(n.mcmc,1,out.1$p.save)
y.pred=rep(0,n.mcmc)
y.pred[z.pred==1]=rpois(sum(z.pred),out.1$lam.save)

plot(table(y)/length(y),type="h")
lines(table(y.pred)/n.mcmc,type="h",lty=1,lwd=4,col=rgb(1,0,0,.5))

###
###  Fit Homogeneous ZINB Model
###

source("zinb.mcmc.R")
n.mcmc=10000
out.2=zinb.mcmc(y,n.mcmc)

layout(matrix(1:3,3,1))
plot(out.2$lam.save,type="l")
plot(out.2$p.save,type="l")
plot(out.2$N.save,type="l")

###
### Posterior Predictive Distribution
###

z.pred=rbinom(n.mcmc,1,out.2$p.save)
y.pred=rep(0,n.mcmc)
y.pred[z.pred==1]=rnbinom(sum(z.pred),mu=out.2$lam.save,size=out.2$N.save)

plot(table(y)/length(y),type="h")
lines(table(y.pred)/n.mcmc,type="h",lty=1,lwd=4,col=rgb(1,0,0,.5))



