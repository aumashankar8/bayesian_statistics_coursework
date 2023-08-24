# set working directory
setwd("~/UT Austin/Spring 2023/Bayesian Stats/hw2/hw2")

# source MCMC for Gibbs sampling
source("mcmc_randomWalk.R")

# read auto data
library(readr)
auto <- read_csv("auto.csv")
#View(auto)

# Setup y and hyperparameters
y <- auto$mpg
mu0 <- mean(y)
sigma20 <- var(y)

# Hyper parameter for log(sigma)
mu_log <- 0
sig_log <- 2.25


# MCMC (Gibbs) for Normal-Normal-Normal(Log)
mcmc.out <- mcmc.randomWalk(y, mu0, sigma20, mu_log, sig_log, sigma2.start = 25, sig.tune = 1, n.mcmc = 50000)

# trace plot
plot(mcmc.out$mu.save, type='l', xlab='k', ylab=bquote(mu), 
     main =bquote("Trace Plot for " ~ mu ))

plot(mcmc.out$log_sigma.save, type='l', xlab='k', ylab=bquote(log(sigma)),
     main =bquote("Trace Plot for " ~ log(sigma)))

#Bayesian Regression Models


#Design Matrices
n <- length(y)

L <- 3
X.list=vector("list",L)
X.list[[1]]=model.matrix(mpg ~ cylinders + model.year, data = auto)
X.list[[2]]=model.matrix(mpg ~ cylinders , data = auto)
X.list[[3]]=model.matrix(mpg ~ model.year, data = auto)

source("regressionDIC.R")
beta.mean.list <- vector("list",L)
beta.var.list  <- vector("list",L)

#list of means and vars for Betas
for(l in 1:L){
  beta.mean.list[[l]]=rep(0,dim(X.list[[l]])[2])
  beta.var.list[[l]]=1000
}

mcmc.out.list=vector("list",L)
DIC.vec=rep(0,L)

for(l in 1:L){
  mcmc.out.list[[l]] <- norm.DIC(y, X.list[[l]], mcmc.out$sigma2.save, beta.mean.list[[l]], beta.var.list[[l]],50000)
  DIC.vec[l]=mcmc.out.list[[l]]$DIC
}

DIC.vec 
#Model 1: 2641.365  #BEST B0: -16.6998405  B1=-3.0051047   B2=0.7448058
#Model 2: 3031.059 
#Model 3: 3980.787
plot(1:L,DIC.vec,type="o",lwd=2,ylab="DIC",xlab="Model")

#Prediction cylinders = 8, model.year = 77
mpg_prediction <- -16.69 + -3.005*8 + 0.7448*77 #16.6196 

