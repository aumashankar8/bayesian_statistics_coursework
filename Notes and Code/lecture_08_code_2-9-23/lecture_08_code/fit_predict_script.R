####
####  Simulate Data 
####

set.seed(2022)
n=30
X=matrix(1,n,2)
X[,2]=rnorm(n)
beta=c(0,1)
s2=.5^2
y=rnorm(n,X%*%beta,sqrt(s2))

plot(X[,2],y)

n.pred=100
X.pred=matrix(1,n.pred,2)
X.pred[,2]=seq(-6,6,,n.pred)

####
####  Fit Model
####

source("norm.reg.2.mcmc.R")
n.mcmc=10000
mcmc.out=norm.reg.2.mcmc(y,X,X.pred,rep(0,2),100,n.mcmc)

####
####  Compare Fitted vs. Predictive Distn  
####

pred.mn=apply(mcmc.out$y.pred.save,1,mean)
pred.CI=t(apply(mcmc.out$y.pred.save,1,quantile,c(0.025,0.975)))
fit.mn=apply(mcmc.out$y.fit.save,1,mean)
fit.CI=t(apply(mcmc.out$y.fit.save,1,quantile,c(0.025,0.975)))

layout(matrix(c(1,1,1,1,2,3),2,3))
matplot(X.pred[,2],pred.CI,type="l",col=2,lwd=2,lty=1,xlab="x",ylab="y")
matplot(X.pred[,2],fit.CI,type="l",col=3,lwd=2,lty=1,add=TRUE)
points(y,X[,2],lwd=2)
legend("topleft",lwd=4,col=2:3,lty=1,legend=c("Prediction CI","Fitted CI"))
plot(pred.CI[,2]-pred.CI[,1],type="l",ylab="pred interval width",xlab="x",col=2)
plot(fit.CI[,2]-fit.CI[,1],type="l",ylab="fit interval width",xlab="x",col=3)


