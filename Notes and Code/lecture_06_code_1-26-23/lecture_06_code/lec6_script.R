####
####  Simulate Data
####

set.seed(2022)
beta.0.true=-2
beta.1.true=1
s2.true=.1

n=30
x=runif(n)
y=rnorm(n,beta.0.true+beta.1.true*x,sqrt(s2.true))

plot(x,y,pch=16,col=rgb(0,0,0,.5))
abline(a=beta.0.true,b=beta.1.true,col=rgb(0,1,0,.5),lwd=2)
legend("topleft",lwd=2,col=rgb(0,1,0,.5),legend=c("Truth"))

####
####  Fit SLR Model 
####

source("slr.mcmc.R")
n.mcmc=50000
out=slr.mcmc(y,x,n.mcmc)

####
####  Trace Plots 
####

layout(matrix(1:3,3,1))
plot(out$beta.0.save,type="l")
plot(out$beta.1.save,type="l")
plot(out$s2.save,type="l")

####
####  Inference 
####

dinvgamma <- function(x,q,r){
  x^(-q-1) * exp(-1/r/x) / (r^q) / gamma(q)
}

layout(matrix(1:3,1,3))
hist(out$beta.0.save,breaks=40,xlab=bquote(beta[0]),main="",prob=TRUE)
curve(dnorm(x,out$mu.0,sqrt(out$s2.0)),col=rgb(1,0,0,.5),lwd=2,add=TRUE)
abline(v=beta.0.true,col=rgb(0,1,0,.4),lwd=2)
hist(out$beta.1.save,breaks=40,xlab=bquote(beta[1]),main="",prob=TRUE)
curve(dnorm(x,out$mu.1,sqrt(out$s2.1)),col=rgb(1,0,0,.5),lwd=2,add=TRUE)
abline(v=beta.1.true,col=rgb(0,1,0,.4),lwd=2)
hist(out$s2.save,breaks=40,xlab=bquote(sigma^2),main="",prob=TRUE)
curve(dinvgamma(x,out$q,out$r),col=rgb(1,0,0,.5),lwd=2,add=TRUE)
abline(v=s2.true,col=rgb(0,1,0,.4),lwd=2)


quantile(out$beta.0.save,c(.025,.975))
quantile(out$beta.1.save,c(.025,.975))

####
####  Fitted Values 
####

n.fitted=100
fitted.mat=matrix(0,n.mcmc,n.fitted)
x.fitted=seq(0,1,,n.fitted)
for(k in 1:n.mcmc){
  fitted.mat[k,]=out$beta.0.save[k]+out$beta.1.save[k]*x.fitted
}

fitted.mn=apply(fitted.mat,2,mean)
fitted.l=apply(fitted.mat,2,quantile,.025)
fitted.u=apply(fitted.mat,2,quantile,.975)

plot(x,y,pch=16,col=rgb(0,0,0,.5))
abline(a=beta.0.true,b=beta.1.true,col=rgb(0,1,0,.5),lwd=2)
polygon(c(x.fitted,rev(x.fitted)),c(fitted.u,rev(fitted.l)),col=rgb(0,0,0,.2),border=NA)
lines(x.fitted,fitted.mn,type="l",lwd=2,col=rgb(0,0,0,.8))
legend("topleft",lwd=2,col=c(rgb(0,1,0,.5),rgb(0,0,0,.8)),legend=c("Truth","Fitted"))



