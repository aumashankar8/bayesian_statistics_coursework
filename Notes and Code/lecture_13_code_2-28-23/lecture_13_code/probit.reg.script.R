####
####  Read in Binary Willow Tit Data (Royle and Dorazio (2008), page 88)
####

wt.df=read.table("wt.1.binary.txt",header=TRUE)

####
####  Create Variables and Scale Covariates 
####

y=wt.df$y
n=length(y)
X1=cbind(rep(1,n),as.matrix(scale(wt.df[,-1])))
X2=cbind(X1,X1[,2]^2)  # include quadratic effect on elevation

####
####  Fit Model with no Quadratic Effect 
####

source("probit.reg.mcmc.R")
wt.1.out=probit.reg.mcmc(y,X1,rep(0,dim(X1)[2]),100,10000)

matplot(t(wt.1.out$beta.save),type="l",lty=1)

boxplot(data.frame(pnorm(t(wt.1.out$z.save))),col=8,outline=FALSE)
lines(y,col=3,lwd=2)

####
####  Fit Model with Quadratic Effect
####

source("probit.reg.mcmc.R")
wt.2.out=probit.reg.mcmc(y,X2,rep(0,dim(X2)[2]),100,10000)

matplot(t(wt.2.out$beta.save),type="l",lty=1)

boxplot(data.frame(pnorm(t(wt.2.out$z.save))),col=8,outline=FALSE)
lines(y,col=3,lwd=2)



