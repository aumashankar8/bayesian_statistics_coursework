####
####  R Script to simulate data and fit normal-normal-inverse-gamma model 
####

set.seed(2021)
mu.true=3
s2.true=2^2
y=rnorm(10,mu.true,sqrt(s2.true))

mu.0=0
s2.0=100
source("norm.IG.mcmc.R")
n.mcmc=10000
n.burn=round(.1*n.mcmc)
out=norm.IG.mcmc(y,mu.0,s2.0,mean(y),n.mcmc)

####
####  Trace Plots
####

layout(matrix(1:2,2,1))
plot(out$mu.save,type="l",main="",ylab=bquote(mu),xlab="MCMC Iteration")
abline(h=mu.true,lwd=2,col=rgb(0,1,0,.5))
plot(out$s2.save,type="l",main="",ylab=bquote(sigma^2),xlab="MCMC Iteration")
abline(h=s2.true,lwd=2,col=rgb(0,1,0,.5))

####
####  Posterior Means and CIs 
####

mean(out$mu.save[-(1:n.burn)])
quantile(out$mu.save[-(1:n.burn)],c(0.025,0.975))

mean(out$s2.save[-(1:n.burn)])
quantile(out$s2.save[-(1:n.burn)],c(0.025,0.975))

####
####  Posterior Histograms 
####

dinvgamma <- function(x,q,r){
  x^(-q-1) * exp(-1/r/x) / (r^q) / gamma(q)
}

layout(matrix(1:2,1,2))
plot(density(out$mu.save[n.burn:n.mcmc]),lwd=2,main="Posterior and Prior: mu")
curve(dnorm(x,mu.0,sqrt(s2.0)),col=2,lwd=2,add=TRUE)
abline(v=mean(y),col=4)
legend("topright",col=c(1,2,4),lwd=2,legend=c("Posterior","Prior","Sample"))
plot(density(out$s2.save[n.burn:n.mcmc]),lwd=2,main="Posterior and Prior: s2")
curve(dinvgamma(x,out$q,out$r),col=2,lwd=2,add=TRUE)
abline(v=var(y),col=4)
legend("topright",col=c(1,2,4),lwd=2,legend=c("Posterior","Prior","Sample"))

