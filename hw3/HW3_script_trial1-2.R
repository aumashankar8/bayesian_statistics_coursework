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
source("hw3-3.betareg.mcmc.R")
n.mcmc=100000

mu.beta=0
s2.beta=1
beta.tune=0.1
gamma1.tau=1
gamma2.tau=4

out.mcmc=hw3.betareg.mcmc(y,X,mu.beta,s2.beta,beta.tune,gamma1.tau,gamma2.tau,n.mcmc)
out.mcmc$p.value #0.99412

#Trace plots
layout(matrix(1:3,3,1))
plot(out.mcmc$beta.save[1,],type="l",xlab= "k", ylab=bquote(b0), 
     main =bquote("Trace Plot for " ~ b0 ))
plot(out.mcmc$beta.save[2,],type="l",xlab= "k", ylab=bquote(b1), 
     main =bquote("Trace Plot for " ~ b1 ))
plot(out.mcmc$beta.save[3,],type="l",xlab= "k", ylab=bquote(b2), 
     main =bquote("Trace Plot for " ~ b2 ))
layout(matrix(1:2,2,1))
plot(out.mcmc$beta.save[4,],type="l",xlab= "k", ylab=bquote(b3), 
     main =bquote("Trace Plot for " ~ b3 ))
plot(out.mcmc$tau.save,type="l",xlab= "k", ylab=bquote(tau), 
     main =bquote("Trace Plot for " ~ tau ))

#Predictions
yprediction=apply(out.mcmc$ypred.save[,-(1:n.burn)],1,mean)
