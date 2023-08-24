####
####  Simulate Data
####

n=10
p.true=.7
set.seed(2022)
y=rbinom(n,1,p.true)
y

####
####  Compile MCMC Algorithm
####

source("binom.beta.mcmc.R")

####
####  Fit Model to Data  
####

alpha=1
beta=1
n.mcmc=10000
binom.beta.out=binom.beta.mcmc(y,alpha,beta,n.mcmc)

####
####  View Trace Plot and M-H acceptance rate 
####

plot(binom.beta.out$p.save,type="l",ylab="p",main=paste("Acceptance: ",binom.beta.out$mh.p/binom.beta.out$n.mcmc))
abline(h=p.true,col=rgb(0,1,0,.5),lwd=3)

####
####  Compute Posterior Summary Statistics
####

mean(binom.beta.out$p.save[-(1:500)])
quantile(binom.beta.out$p.save[-(1:500)],c(0.025,0.975))

####
####  Histogram of MCMC Output
####

hist(binom.beta.out$p.save,breaks=50,xlim=c(0,1),prob=TRUE,col=8,main=expression(paste("Posterior for ",p)),xlab="p")
curve(dbeta(x,sum(y)+alpha,sum(1-y)+beta),lwd=2,add=TRUE)
abline(v=p.true,col=rgb(0,1,0,.5),lwd=3)
abline(v=mean(binom.beta.out$p.save[-(1:500)]),col=rgb(1,0,0,.5),lwd=3)
abline(v=quantile(binom.beta.out$p.save[-(1:500)],c(0.025,0.975)),col=rgb(1,0,0,.5),lwd=3,lty=2)
legend("topleft",lty=c(1,1,2),lwd=3,col=c(3,2,2),legend=c("truth","post mean","CI"))


