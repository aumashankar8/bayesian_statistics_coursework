###
###  MC approx unif mean
###

K=100000
theta=runif(K)

mean(theta)

hist(theta,breaks=50,prob=TRUE,col=8,xlab=bquote(theta),main="",ylim=c(0,1.25))
curve(dunif(x),col=3,lwd=2,add=TRUE)
legend("topright",fill=c(8,NA),col=c(NA,3),border=c(1,NA),lwd=c(NA,2),x.intersp=c(-.5,1),bty="n",legend=c("Monte Carlo","Exact"))

###
###  mean and distribution of transformation
###

mean(theta^2)

hist(theta^2,breaks=50,prob=TRUE,xlab=bquote(theta^2),main="")
curve(1/2/sqrt(x),col=3,lwd=2,add=TRUE)
legend("topright",fill=c(8,NA),col=c(NA,3),border=c(1,NA),lwd=c(NA,2),x.intersp=c(-.5,1),bty="n",legend=c("Monte Carlo","Exact"))

###
###  joint vs marginal distribution
###

K=100000
theta.1=runif(K,0,1)
theta.2=runif(K,0,theta.1)
cor(theta.1,theta.2)

layout(matrix(1:4,2,2))
hist(theta.2,breaks=50,prob=TRUE,col=8,xlab=bquote(theta[2]),main="")
curve(-log(x),col=3,lwd=2,add=TRUE)
legend("topright",fill=c(8,NA),col=c(NA,3),border=c(1,NA),lwd=c(NA,2),x.intersp=c(-.5,1),bty="n",legend=c("Monte Carlo","Exact"))
plot.new()
plot(theta.1,theta.2,pch=16,cex=.2,col=rgb(0,0,0,.25),xlim=c(0,1),ylim=c(0,1),xlab=bquote(theta[1]),ylab=bquote(theta[2]))
hist(theta.1,breaks=50,prob=TRUE,col=8,xlab=bquote(theta[1]),main="",ylim=c(0,1.25))
curve(dunif(x),col=3,lwd=2,add=TRUE)
legend("topright",fill=c(8,NA),col=c(NA,3),border=c(1,NA),lwd=c(NA,2),x.intersp=c(-.5,1),bty="n",legend=c("Monte Carlo","Exact"))

mean(theta.1)
mean(theta.2)


