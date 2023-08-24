####
####  Read in Binary Willow Tit Data (Royle and Dorazio (2008), page 88)
####

wt.df=read.table("wt.1.binary.txt",header=TRUE)

####
####  Create Variables and Scale Covariates 
####

y=wt.df$y
n=length(y)
X1=cbind(rep(1,n),as.matrix(scale(wt.df[,-1]))) # intercept, elevation, forest
X2=cbind(X1,X1[,2]^2) # include quadratic effect on elevation 

####
####  Fit Model with no Quadratic Effect 
####

source("logit.reg.mcmc.R")
wt.1.out=logit.reg.mcmc(y,X1,rep(0,dim(X1)[2]),1.5,.1,10000)

matplot(t(wt.1.out$betasave),type="l",lty=1)

####
####  Fit Model with Quadratic Effect 
####

source("logit.reg.mcmc.R")
wt.2.out=logit.reg.mcmc(y,X2,rep(0,dim(X2)[2]),1.5,.1,10000)

matplot(t(wt.2.out$betasave),type="l",lty=1)



