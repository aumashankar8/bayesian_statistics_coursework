####
####  Read in Data
####

med.df=read.csv("MedicalData.csv",header=TRUE)
head(med.df)

####
####  View Data   
####

y=med.df$BP
n=length(y)
X=matrix(1,n,3)
X[,2]=med.df$BMI
X[,3]=med.df$Glucose

pairs(cbind(y,X[,-1])) 
cor(cbind(y,X[,-1]))

####
####  Fit Model 
####

source("norm.reg.mcmc.R")
n.mcmc=10000
beta.mn=rep(0,3)
beta.var=1000
mcmc.out=norm.reg.mcmc(y,X,beta.mn,beta.var,n.mcmc,TRUE)

####
####  Make Trace Plots 
####

layout(matrix(1:4,2,2))
plot(mcmc.out$beta.save[1,],type="l",lty=1,ylab=bquote(beta[0]))
plot(mcmc.out$beta.save[2,],type="l",lty=1,ylab=bquote(beta[1]))
plot(mcmc.out$beta.save[3,],type="l",lty=1,ylab=bquote(beta[2]))
plot(mcmc.out$s2.save,type="l",ylab=bquote(sigma^2))

####
####  Make Posterior Histograms 
####

dIG <- function(x,q,r){
  x^(-q-1) * exp(-1/r/x) / (r^q) / gamma(q)
}

layout(matrix(1:4,2,2))
hist(mcmc.out$beta.save[1,],prob=TRUE,breaks=40,xlab=bquote(beta[0]),main="")
curve(dnorm(x,beta.mn[1],sqrt(beta.var)),lwd=2,add=TRUE,col=rgb(1,0,0,.5))
hist(mcmc.out$beta.save[2,],prob=TRUE,breaks=40,xlab=bquote(beta[1]),main="")
curve(dnorm(x,beta.mn[2],sqrt(beta.var)),lwd=2,add=TRUE,col=rgb(1,0,0,.5))
hist(mcmc.out$beta.save[3,],prob=TRUE,breaks=40,xlab=bquote(beta[2]),main="")
curve(dnorm(x,beta.mn[3],sqrt(beta.var)),lwd=2,add=TRUE,col=rgb(1,0,0,.5))
hist(mcmc.out$s2.save,prob=TRUE,breaks=40,xlab=bquote(sigma^2),main="")
curve(dIG(x,mcmc.out$q,mcmc.out$r),lwd=2,add=TRUE,col=rgb(1,0,0,.5))

