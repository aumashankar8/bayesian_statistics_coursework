#Import data
library(betareg)
data("GasolineYield", package="betareg")
crude.df <- GasolineYield
crude.df

#Standardize covariates or dependent variables
crude.df$gravity=((crude.df$gravity-mean(crude.df$gravity))^2)/sqrt(var(crude.df$gravity))
crude.df$pressure=((crude.df$pressure-mean(crude.df$pressure))^2)/sqrt(var(crude.df$pressure))
crude.df$temp=((crude.df$temp-mean(crude.df$temp))^2)/sqrt(var(crude.df$temp))

crude.df

#Create data vector and design matrix
y=crude.df$yield
X=model.matrix(yield~gravity+pressure+temp,data=crude.df)

N=dim(crude.df)[1]
p=dim(X)[2]

#Fit Beta regression model
source("hw3.betareg.mcmc.R")
n.mcmc=10000

mu.beta=rep(0,p)
s2.beta=100^2
beta.tune=1
gamma1.tau=0.1
gamma2.tau=0.2

out.mcmc=hw3.betareg.mcmc(y,X,mu.beta,s2.beta,beta.tune,gamma1.tau,gamma2.tau,n.mcmc)

#Trace plots
layout(matrix(1:4,4,1))
plot(out.mcmc$beta.save[1,],type="l",lty=1)
plot(out.mcmc$beta.save[2,],type="l",lty=1)
plot(out.mcmc$beta.save[3,],type="l",lty=1)
plot(out.mcmc$tau.save[1,],type="l",lty=1)
