#### 
####  R code for creating mixture figures 
#### 

curve(dnorm(x,-2,1),col=2,lwd=3,xlim=c(-5,7),ylab="density",xlab="y")
curve(dnorm(x,1,2),col=3,lwd=3,add=TRUE)

p=.75
curve(dnorm(x,-2,1),col=rgb(1,0,0,.5),lwd=3,xlim=c(-5,7),ylab="density",xlab="y",main="p = 0.75")
curve(dnorm(x,1,2),col=rgb(0,1,0,.5),lwd=3,add=TRUE)
curve(p*dnorm(x,-2,1)+(1-p)*dnorm(x,1,2),col=8,lwd=5,add=TRUE)

p=.5
curve(dnorm(x,-2,1),col=rgb(1,0,0,.5),lwd=3,xlim=c(-5,7),ylab="density",xlab="y",main="p = 0.5")
curve(dnorm(x,1,2),col=rgb(0,1,0,.5),lwd=3,add=TRUE)
curve(p*dnorm(x,-2,1)+(1-p)*dnorm(x,1,2),col=8,lwd=5,add=TRUE)

p=.25
curve(dnorm(x,-2,1),col=2,lwd=3,xlim=c(-5,7),ylab="density",xlab="y",main="p = 0.25")
curve(dnorm(x,1,2),col=3,lwd=3,add=TRUE)
curve(p*dnorm(x,-2,1)+(1-p)*dnorm(x,1,2),col=8,lwd=5,add=TRUE)

###
###  Read in Old Faithful Data
###

faithful.df=read.csv("faithful.csv",header=TRUE)
head(faithful.df)

y=faithful.df$Eruption.wait..mins.
n=length(y)
hist(y,col=rgb(0,0,0,.4),main="",xlab="y (mins)")

###
###  Fit Mixture Model to Data 
###

source("mix.mcmc.R")
n.mcmc=10000
out=mix.mcmc(y,n.mcmc)

layout(matrix(1:3,3,1))
matplot(out$mu.save,type="l",lty=1)
matplot(out$s2.save,type="l",lty=1)
plot(out$p.save,type="l")

###
###  Graph Results including PPD
###

z.pred=rbinom(n.mcmc,1,out$p.save)
y.pred=z.pred*rnorm(n.mcmc,out$mu.save[,1],sqrt(out$s2.save[,1]))+(1-z.pred)*rnorm(n.mcmc,out$mu.save[,2],sqrt(out$s2.save[,2]))

d.1=density(out$mu.save[,1])
d.2=density(out$mu.save[,2])
d.pred=density(y.pred)

hist(y,col=rgb(0,0,0,.4),main="",xlab="y (mins)",prob=TRUE,ylim=c(0,max(d.1$y,d.2$y)),breaks=20)
lines(d.1,col=rgb(1,0,0,.5),lwd=2)
lines(d.2,col=rgb(0,1,0,.5),lwd=2)
lines(d.pred,col=rgb(0,0,1,.5),lwd=2)

