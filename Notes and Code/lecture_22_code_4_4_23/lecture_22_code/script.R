###
###  Read in Data
###

dat.df=read.csv("oxboys.csv",header=TRUE)

n=max(dat.df$Subject)
J=sum(dat.df$Subject==1)

###
###  Display Scatterplot 
###

plot(dat.df$age,dat.df$height,type="p",ylab="height",xlab="age")

###
###  Fit Hierarchical Model w/ Random Coefficients 
###

source("hier.mcmc.R")
n.mcmc=10000
out=hier.mcmc(dat.df$height,dat.df$age,dat.df$Subject,n.mcmc)

###
###  Posterior CIs and Means 
###

all.0=rbind(out$beta.save[1,,],out$mu.save[1,])
all.1=rbind(out$beta.save[2,,],out$mu.save[2,])

layout(matrix(1:2,2,1))
plot(1:(n+1),rep(0,n+1),ylim=range(all.0),ylab=bquote(beta["0,i"]),xlab="individual",xaxt="n",type="n")
mtext(bquote(mu[0]),4,.5,cex=1.25)
axis(1,1:(n+1),labels=c(1:n,"pop"),cex.axis=1.1)
segments(1:(n+1),apply(all.0,1,quantile,0.025),1:(n+1),apply(all.0,1,quantile,0.975),lwd=2,lty=c(rep(1,n),6))
points(1:(n+1),apply(all.0,1,mean),pch=20,cex=.5)
abline(h=mean(out$mu.save[1,]),col=8)
plot(1:(n+1),rep(0,n+1),ylim=range(all.1),ylab=bquote(beta["1,i"]),xlab="individual",xaxt="n",type="n")
mtext(bquote(mu[1]),4,.5,cex=1.25)
axis(1,1:(n+1),labels=c(1:n,"pop"),cex.axis=1.1)
segments(1:(n+1),apply(all.1,1,quantile,0.025),1:(n+1),apply(all.1,1,quantile,0.975),lwd=2,lty=c(rep(1,n),6))
points(1:(n+1),apply(all.1,1,mean),pch=20,cex=.5)
abline(h=mean(out$mu.save[2,]),col=8)

###
###  Posterior Mean Regression Lines 
###

plot(dat.df$age,dat.df$height,type="p",ylab="height",xlab="age")
for(i in 1:n){
  abline(mean(out$beta.save[1,i,]),mean(out$beta.save[2,i,]),col=rgb(0,0,0,.5))
}

###
###  Posterior histogram for rho (correlation) 
###

hist(out$rho.save,breaks=30,xlim=c(-1,1))
mean(out$rho.save>0)



