####
####  First Try with Simulations
####

source("ar1.sim.mcmc.R")
tmp.out=ar1.sim.mcmc(100,.9,10000)
#tmp.out=ar1.sim.mcmc(100,.99,10000)

####
####  Load Scaup Data 
####

st30.df=read.table("st30.dat",header=TRUE)

####
####  Plot Scaup Data 
####

layout(matrix(1:2,2,1))
plot(st30.df$year,st30.df$y,type="l",lty=1,lwd=2)
plot(st30.df$year,log(st30.df$y),type="l",lty=1,lwd=2)

####
####  Fit AR(1) Model to Stratum 30 log Scaup Data 
####

y=log(st30.df$y)
source("ar1.mcmc.R")
n.mcmc=10000
out=ar1.mcmc(y,n.mcmc)

####
####  One step ahead forecast 
####

T=length(y)
y.pred.save=rnorm(n.mcmc,out$alpha.save*(y[T]-mean(y)),sqrt(out$s2.save))+mean(y)

y.pred.mean=mean(y.pred.save)
hist(y.pred.save,breaks=40,main="")

plot(c(st30.df$year,max(st30.df$year)+1),c(y,y.pred.mean),type="l",lty=1,lwd=2)
points(max(st30.df$year)+1,y.pred.mean,col=2,pch=16)

