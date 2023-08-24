#Import data
library(betareg)
data("GasolineYield", package="betareg")
crude.df <- GasolineYield
crude.df

#Standardize covariates or dependent variables
crude.df$gravity=scale(crude.df$gravity)
crude.df$pressure=scale(crude.df$pressure)
crude.df$temp=scale(crude.df$temp)

crude.df

#Create data vector and design matrix
y=crude.df$yield
X=model.matrix(yield~gravity+pressure+temp,data=crude.df)

N=dim(crude.df)[1]
p=dim(X)[2]

#Fit Beta regression model
source("hw3-4.betareg.mcmc.R")
n.mcmc=100000

mu.beta=0
s2.beta=0.075
beta.tune=0.15
gamma1.tau=1
gamma2.tau=4

out.mcmc=hw3.betareg.mcmc(y,X,mu.beta,s2.beta,beta.tune,gamma1.tau,gamma2.tau,n.mcmc)

#Trace plots
layout(matrix(1:2,2,1))
plot(out.mcmc$beta.save[1,],main="intercept",ylab=bquote(beta[0]),xlab="K",type="l",lty=1)
plot(out.mcmc$beta.save[2,],main="gravity",ylab=bquote(beta[1]),xlab="K",type="l",lty=1)

layout(matrix(1:2,2,1))
plot(out.mcmc$beta.save[3,],main="pressure",ylab=bquote(beta[2]),xlab="K",type="l",lty=1)
plot(out.mcmc$beta.save[4,],main="temp",ylab=bquote(beta[3]),xlab="K",type="l",lty=1)

plot(out.mcmc$tau.save,ylab=bquote(tau),xlab="K",type="l",lty=1)
