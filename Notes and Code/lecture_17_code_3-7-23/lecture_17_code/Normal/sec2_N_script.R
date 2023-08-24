###
###  Read in T_max Data
###
#
#  Troia, M.J. and Giam, X., 2019. Extreme heat events and the vulnerability of endemic montane fishes to climate change. Ecography, 42(11), pp.1913-1925.
#
setwd("~/UT Austin/Spring 2023/Bayesian Stats/Notes and Code/lecture_17_code/lecture_17_code/Normal")
tmax.df=read.csv("laboratory_CTM.csv")
head(tmax.df)

n=dim(tmax.df)[1]
y=c(tmax.df$LRR_sustained)

p=5
X=matrix(1,n,p)
X[,1:4]=model.matrix(~0+tmax.df$species) #Dummy variable for each species. 1 or a 0. 4 different speciies
                                         # Traditional way is to use 1 intercept and then 3 other variables that are relative to that BASE species.
X[,5]=tmax.df$acclim_temp

pairs(cbind(y,X[,-1]))
summary(lm(y~0+X))  # check least squares fit 

###
###  Specify Linear Regression Model in JAGS
###

library(rjags)

m.jags <-"
  model{
    for(i in 1:n){
      y[i] ~ dnorm(mu[i], tau)  # note: tau is precision=1/variance!  
      mu[i] = X[i,1:p]%*%beta
    }
    beta ~ dmnorm(mu.beta,Tau.beta)
    tau <- pow(sig, -2)
    sig ~ dunif(0,10)
}
"
#pow()=power or exponent
mod<-textConnection(m.jags)

###
###  Fit Linear Regression Model w/ JAGS 
###

mu.beta=rep(0,p)
Sig.beta=10000*diag(p)
Tau.beta=solve(Sig.beta)
m.out<-jags.model(mod,data=list('y'=y,'n'=n,'X'=X,'p'=p,'mu.beta'=mu.beta,'Tau.beta'=Tau.beta),n.chains=1,n.adapt=0) # build model and algorithm

n.mcmc=50000
n.burn=round(.2*n.mcmc)
update(m.out,n.burn) # perform burn-in
m.samples=jags.samples(m.out,c('beta','sig'),n.mcmc) # fit model post-burn-in

layout(matrix(1:2,2,1))
matplot(t(m.samples$beta[,,1]),type="l",lty=1,ylab=bquote(beta))
plot(m.samples$sig[1,,1],type="l",lty=1,ylab=bquote(sigma))

###
###  Linear Regression Inference 
###

library(vioplot)
layout(matrix(c(1,1,1,2),1,4))
vioplot(data.frame(t(m.samples$beta[1:4,,1])),names=expression(beta[1],beta[2],beta[3],beta[4]))
vioplot(data.frame(m.samples$beta[5,,1]),names=expression(beta[5]))
abline(h=0,col=8)

apply(m.samples$beta[,,1],1,mean) # marginal posterior means for beta 
apply(m.samples$beta[,,1],1,sd) # marginal posterior sd for beta 
apply(m.samples$beta[,,1],1,quantile,c(0.025,.975)) # marginal posterior 95% CI for beta 

mean(m.samples$beta[1,,1]>m.samples$beta[2,,1]) # Post prob: P(beta_1 > beta_2 | y) = 0.99
plot(m.samples$beta[1,,1],m.samples$beta[2,,1],col=rgb(0,0,0,.1),pch=16,cex=.5,asp=TRUE) # don't forget about joint inference! 
abline(0,1,col=2)

